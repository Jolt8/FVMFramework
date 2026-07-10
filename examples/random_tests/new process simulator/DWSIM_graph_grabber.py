import json
import clr
import sys
import os

# 1. Add the DWSIM installation path to sys.path
dwsimpath = r"C:\Users\wille\AppData\Local\DWSIM"
sys.path.append(dwsimpath)

# 2. Add references to DWSIM assemblies
clr.AddReference(os.path.join(dwsimpath, "DWSIM.Automation.dll"))
clr.AddReference(os.path.join(dwsimpath, "DWSIM.Interfaces.dll"))

# 3. Now Python will recognize the namespaces
from DWSIM.Interfaces.Enums import PropertyType
from DWSIM.Automation import Automation3

# 4. Initialize automation and load flowsheet
interf = Automation3()

# Get file path from command line arguments or use default
if len(sys.argv) > 1:
    filepath = sys.argv[1]
else:
    # Use absolute path so it can be run from anywhere
    filepath = os.path.join(os.path.dirname(os.path.abspath(__file__)), "small_recycle_test.dwxmz")

print(f"Loading flowsheet from: {filepath}")
Flowsheet = interf.LoadFlowsheet(filepath)

if Flowsheet is None:
    print("Failed to load flowsheet.")
    sys.exit(1)
edges = []

#print(Flowsheet.SimulationObjects) #this has almost everything we need except connectivity
#print(dir(Flowsheet.SimulationObjects))
#print(Flowsheet.SimulationObjects.Values[1].GetConnectionPortsInfo)
#print(Flowsheet.SimulationObjects.Values.GetConnectionPortsList)

#print(dir(Flowsheet.SimulationObjects.Values.GraphicObject))

for obj in Flowsheet.SimulationObjects.Values:
    # Safely inspect the GraphicObject layer
    if hasattr(obj, "GraphicObject") and obj.GraphicObject is not None:
        g_obj = obj.GraphicObject

        #print(dir(g_obj))
        #print("info", obj.GetConnectionPortsInfo())
        #print("list", obj.GetConnectionPortsList())

        for connection in obj.GetConnectionPortsInfo():
            if connection.IsConnected == True:
                print("connection", connection)
                print("connected object", connection.ConnectedObject)
                print("name", connection.Name)
                print("tag", connection.ConnectedObject.GraphicObject.Tag)
                print("category", connection.ConnectedObject.ObjectClass)
                print("type", connection.ConnectedObject.GetType().Name)
                #print("properties", connection.ConnectedObject.GetProperties(PropertyType.ALL))
                #print("properties 2", connection.ConnectedObject.GetProperties2())
                #print("properties", connection.ConnectedObject.ExtraProperties)

                if str(connection.ConnectedObject.GetType().Name) != "MaterialStream":
                    sim_obj_interface = connection.ConnectedObject
                    sim_obj = sim_obj_interface.GetAsObject()
                    print(f"--- Properties for {sim_obj.GraphicObject.Tag} ---")
                    # 1. Get the list of available property names
                    prop_names = sim_obj.GetProperties(PropertyType.ALL)
                    # 2. Loop through them and get their values
                    for prop_name in prop_names:
                        # GetPropertyValue uses a C# "ref" parameter. 
                        # In Python, this returns a tuple: (success_boolean, value_string)
                        prop_value = sim_obj.GetPropertyValue(prop_name)
                                
                        # Optionally, you can also get the unit of measure!
                        prop_unit = sim_obj.GetPropertyUnit(prop_name)
                        print(f"{prop_name}: {prop_value} {prop_unit}")

                    if hasattr(sim_obj, "CalcMode"):
                        mode_enum = sim_obj.CalcMode
                        mode_name = str(mode_enum) # Converts the .NET Enum to a readable string
                        print("CalcMode:", mode_name)

        """
        for conn in g_obj.GetInputConnectors():
            stream = conn.getConnectedObject()
            print("Connected Inlet Stream: " + stream.Name)

        for conn in g_obj.GetOutputConnectors():
            stream = conn.getConnectedObject()
            print("Connected Outlet Stream: " + stream.Name)
        """
        
        # Check if this object connects to input/output streams
        # DWSIM's graphic object interface maps input and output connection ports
        if hasattr(g_obj, "InputConnectors") and hasattr(g_obj, "OutputConnectors"):
            # Connect incoming streams to this Unit Op
            for input_stream_id in g_obj.InputConnectors:
                if input_stream_id:
                    edges.append({"source": str(input_stream_id), "target": str(g_obj.Tag)})
                    
            # Connect this Unit Op to outgoing streams
            for output_stream_id in g_obj.OutputConnectors:
                if output_stream_id:
                    edges.append({"source": str(g_obj.Tag), "target": str(output_stream_id)})
"""
# Format and save JSON
json_data = json.dumps(edges, indent=4)
#print(json_data)

with open("C:/temp/dwsim_graph.json", "w") as f:
    f.write(json_data)
"""




#check: https://dwsim.org/api_help/html/T_DWSIM_SharedClasses_ConnectionPortInfo.htm