"""
extract_dwsim_flowsheet.py

Opens a DWSIM .dwxmz flowsheet via the Automation3 API and extracts:
  - Property package name
  - Compound data (molecular weight, ideal-gas Cp constants, critical properties)
  - Material streams (mass flow, temperature, pressure, mass fractions, vapor fraction)
  - Energy streams (power)
  - Unit operations (all DWSIM properties + CalcMode)
  - Directed edges (connectivity graph)

Writes everything to a single JSON file.

Usage:
    python extract_dwsim_flowsheet.py <path_to_dwxmz> [output_json_path]
"""

import json
import sys
import os

import clr

# ── DWSIM bootstrap ────────────────────────────────────────────────────────────
dwsimpath = r"C:\Users\wille\AppData\Local\DWSIM"
sys.path.append(dwsimpath)

clr.AddReference(os.path.join(dwsimpath, "DWSIM.Automation.dll"))
clr.AddReference(os.path.join(dwsimpath, "DWSIM.Interfaces.dll"))

from DWSIM.Interfaces.Enums import PropertyType
from DWSIM.Automation import Automation3

# ── Helpers ─────────────────────────────────────────────────────────────────────

STREAM_TYPES = {"MaterialStream", "EnergyStream"}


def safe_float(val, default = None):
    """Try to convert a value to float, returning *default* on failure."""
    try:
        return float(val)
    except (TypeError, ValueError):
        return default


def safe_str(val, default = ""):
    """Convert a .NET object to a Python string safely."""
    try:
        return str(val)
    except Exception:
        return default


# ── Extraction functions ────────────────────────────────────────────────────────

def extract_property_packages(flowsheet):
    """Return a list of property-package display names attached to the flowsheet."""
    names = []
    try:
        for pp in flowsheet.PropertyPackages.Values:
            names.append(safe_str(pp.Tag) or safe_str(pp.Name))
    except Exception:
        pass
    return names


def extract_compounds(flowsheet):
    """
    Pull per-compound constant data from Flowsheet.SelectedCompounds.
    Returns a list of dicts.
    """
    compounds = []
    try:
        for comp in flowsheet.SelectedCompounds.Values:
            entry = {
                "name": safe_str(comp.Name),
                "molecular_weight": safe_float(comp.Molar_Weight),  # g/mol
                "critical_temperature": safe_float(comp.Critical_Temperature),  # K
                "critical_pressure": safe_float(comp.Critical_Pressure),  # Pa
                "acentric_factor": safe_float(comp.Acentric_Factor),
            }

            # Ideal-gas Cp constants (for polynomial / equation evaluation)
            for suffix in ("A", "B", "C", "D", "E"):
                attr = f"Ideal_Gas_Heat_Capacity_Const_{suffix}"
                entry[f"ig_cp_{suffix.lower()}"] = safe_float(getattr(comp, attr, None))

            # Equation type used for the Cp correlation
            entry["ig_cp_equation"] = safe_str(getattr(comp, "IdealgasCpEquation", ""))

            compounds.append(entry)
    except Exception as exc:
        print(f"[WARN] Could not extract compounds: {exc}")

    return compounds


def extract_material_stream(sim_obj):
    """
    Extract the fields we care about from a MaterialStream simulation object:
    mass_flow (kg/s), temperature (K), pressure (Pa), mass_fractions, vapor_fraction.
    """
    tag = safe_str(sim_obj.GraphicObject.Tag)

    # Overall phase properties
    overall = sim_obj.GetPhase("Overall")
    props = overall.Properties

    mass_flow = safe_float(props.massflow)           # kg/s
    temperature = safe_float(props.temperature)      # K
    pressure = safe_float(props.pressure)            # Pa
    vapor_fraction = safe_float(getattr(props, "vaporFraction", None))

    # Mass fractions per compound
    mass_fractions = {}
    try:
        for comp in overall.Compounds.Values:
            mass_fractions[safe_str(comp.Name)] = safe_float(comp.MassFraction, 0.0)
    except Exception:
        pass

    return {
        "tag": tag,
        "mass_flow": mass_flow,
        "temperature": temperature,
        "pressure": pressure,
        "vapor_fraction": vapor_fraction,
        "mass_fractions": mass_fractions,
    }


def extract_energy_stream(sim_obj):
    """
    Extract power from an EnergyStream.
    DWSIM stores energy flow in W under property key PROP_ES_0 (or via EnergyFlow).
    """
    tag = safe_str(sim_obj.GraphicObject.Tag)

    power = None
    try:
        # Try the direct attribute first
        power = safe_float(getattr(sim_obj, "EnergyFlow", None))
    except Exception:
        pass

    if power is None:
        try:
            power = safe_float(sim_obj.GetPropertyValue("PROP_ES_0"))
        except Exception:
            pass

    return {
        "tag": tag,
        "power": power,
    }


def extract_unit_operation(sim_obj):
    """
    Extract all properties from a unit operation plus its CalcMode.
    """
    tag = safe_str(sim_obj.GraphicObject.Tag)
    type_name = safe_str(sim_obj.GetType().Name)

    # CalcMode (enum → string)
    calc_mode = None
    actual_obj = sim_obj
    try:
        actual_obj = sim_obj.GetAsObject()
    except Exception:
        pass

    if hasattr(actual_obj, "CalcMode"):
        calc_mode = safe_str(actual_obj.CalcMode)

    # All available properties
    properties = {}
    try:
        prop_names = actual_obj.GetProperties(PropertyType.ALL)
        for prop_name in prop_names:
            pn = safe_str(prop_name)
            try:
                pval = safe_float(actual_obj.GetPropertyValue(prop_name))
                punit = safe_str(actual_obj.GetPropertyUnit(prop_name))
                properties[pn] = {"value": pval, "unit": punit}
            except Exception:
                pass
    except Exception:
        pass

    return {
        "tag": tag,
        "type_name": type_name,
        "calc_mode": calc_mode,
        "properties": properties,
    }


def extract_edges(flowsheet):
    """
    Walk every SimulationObject and build directed edges via GraphicObject connectors.
    Each edge is {source: tag, target: tag}.
    """
    edges = []

    for obj in flowsheet.SimulationObjects.Values:
        g = getattr(obj, "GraphicObject", None)
        if g is None:
            continue

        obj_tag = safe_str(g.Tag)
        obj_type = safe_str(obj.GetType().Name)

        # For streams: an input connector whose AttachedConnector.AttachedFrom is the
        # upstream unit op, and an output connector whose AttachedConnector.AttachedTo
        # is the downstream unit op.
        # For unit ops: input connectors reference incoming streams, output connectors
        # reference outgoing streams.

        # We adopt the convention: stream nodes sit on edges between unit-op nodes.
        # So the connectivity we record is:
        #   unit_op  -->  stream   (unit op output connector points to stream)
        #   stream   -->  unit_op  (stream output connector points to unit op)

        if obj_type in STREAM_TYPES:
            # Input connector → who feeds this stream?
            try:
                for ic in g.InputConnectors:
                    if ic is not None and ic.IsAttached:
                        src = ic.AttachedConnector.AttachedFrom
                        if src is not None:
                            src_tag = safe_str(src.Tag)
                            if src_tag and src_tag != obj_tag:
                                edges.append({"source": src_tag, "target": obj_tag})
            except Exception:
                pass

            # Output connector → who consumes this stream?
            try:
                for oc in g.OutputConnectors:
                    if oc is not None and oc.IsAttached:
                        dst = oc.AttachedConnector.AttachedTo
                        if dst is not None:
                            dst_tag = safe_str(dst.Tag)
                            if dst_tag and dst_tag != obj_tag:
                                edges.append({"source": obj_tag, "target": dst_tag})
            except Exception:
                pass

    # Deduplicate (same edge may appear twice when traversing from both ends)
    seen = set()
    unique_edges = []
    for e in edges:
        key = (e["source"], e["target"])
        if key not in seen:
            seen.add(key)
            unique_edges.append(e)

    return unique_edges


# ── Main ────────────────────────────────────────────────────────────────────────

def main():
    # ── Parse arguments ─────────────────────────────────────────────────────
    if len(sys.argv) < 2:
        print("Usage: python extract_dwsim_flowsheet.py <path_to_dwxmz> [output_json]")
        sys.exit(1)

    filepath = sys.argv[1]
    output_path = sys.argv[2] if len(sys.argv) > 2 else os.path.splitext(filepath)[0] + "_flowsheet_data.json"

    # ── Load flowsheet ──────────────────────────────────────────────────────
    print(f"Loading flowsheet from: {filepath}")
    interf = Automation3()
    flowsheet = interf.LoadFlowsheet(filepath)

    if flowsheet is None:
        print("ERROR: Failed to load flowsheet.")
        sys.exit(1)

    print("Flowsheet loaded successfully.")

    # ── Property packages ───────────────────────────────────────────────────
    property_packages = extract_property_packages(flowsheet)
    print(f"Property packages: {property_packages}")

    # ── Compounds ───────────────────────────────────────────────────────────
    compounds = extract_compounds(flowsheet)
    print(f"Compounds found: {[c['name'] for c in compounds]}")

    # ── Classify simulation objects ─────────────────────────────────────────
    material_streams = []
    energy_streams = []
    unit_operations = []

    for pre_obj in flowsheet.SimulationObjects.Values:
        obj = pre_obj.GetAsObject()

        type_name = safe_str(obj.GetType().Name)

        if type_name == "MaterialStream":
            material_streams.append(extract_material_stream(obj))
        elif type_name == "EnergyStream":
            energy_streams.append(extract_energy_stream(obj))
        else:
            unit_operations.append(extract_unit_operation(obj))

    print(f"Material streams: {[s['tag'] for s in material_streams]}")
    print(f"Energy streams:   {[s['tag'] for s in energy_streams]}")
    print(f"Unit operations:  {[u['tag'] for u in unit_operations]}")

    # ── Edges ───────────────────────────────────────────────────────────────
    edges = extract_edges(flowsheet)
    print(f"Edges found: {len(edges)}")

    # ── Assemble output ─────────────────────────────────────────────────────
    data = {
        "property_packages": property_packages,
        "compounds": compounds,
        "material_streams": material_streams,
        "energy_streams": energy_streams,
        "unit_operations": unit_operations,
        "edges": edges,
    }

    with open(output_path, "w", encoding = "utf-8") as f:
        json.dump(data, f, indent = 4, default = str)

    print(f"\nFlowsheet data written to: {output_path}")


if __name__ == "__main__":
    main()
