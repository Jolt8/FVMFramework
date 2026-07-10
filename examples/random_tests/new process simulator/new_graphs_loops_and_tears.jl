using Graphs
using MetaGraphsNext
using Makie
using GLMakie
using GraphMakie
using NetworkLayout

abstract type AbstractComponent{CalcMode} end

struct Compressor{CalcMode} <: AbstractComponent{CalcMode}
    id::Symbol
    stream_in::Symbol
    stream_out::Symbol
    energy_in::Symbol
end

get_inlet_material_streams(c::Compressor{CalcMode}) where {CalcMode} = (c.stream_in, )
get_inlet_energy_streams(c::Compressor{CalcMode}) where {CalcMode} = (c.energy_in, )
get_outlet_material_streams(c::Compressor{CalcMode}) where {CalcMode} = (c.stream_out, )
get_outlet_energy_streams(c::Compressor{CalcMode}) where {CalcMode} = ()

struct TwoMixer{CalcMode} <: AbstractComponent{CalcMode}
    id::Symbol
    stream_in_1::Symbol
    stream_in_2::Symbol
    stream_out::Symbol
end

get_inlet_material_streams(c::TwoMixer{CalcMode}) where {CalcMode} = (c.stream_in_1, c.stream_in_2)
get_inlet_energy_streams(c::TwoMixer{CalcMode}) where {CalcMode} = ()
get_outlet_material_streams(c::TwoMixer{CalcMode}) where {CalcMode} = (c.stream_out, )
get_outlet_energy_streams(c::TwoMixer{CalcMode}) where {CalcMode} = ()

struct TwoSplitter{CalcMode} <: AbstractComponent{CalcMode}
    id::Symbol
    stream_in::Symbol
    stream_out_1::Symbol
    stream_out_2::Symbol
end

get_inlet_material_streams(c::TwoSplitter{CalcMode}) where {CalcMode} = (c.stream_in, )
get_inlet_energy_streams(c::TwoSplitter{CalcMode}) where {CalcMode} = ()
get_outlet_material_streams(c::TwoSplitter{CalcMode}) where {CalcMode} = (c.stream_out_1, c.stream_out_2)
get_outlet_energy_streams(c::TwoSplitter{CalcMode}) where {CalcMode} = ()

struct TearPoint
    component_1::Symbol
    component_2::Symbol
end

abstract type AbstractStream end

struct MaterialStream <: AbstractStream 
    id::Symbol
end

struct EnergyStream <: AbstractStream 
    id::Symbol
end

unordered_topo = (
    mix_1 = TwoMixer{:AveragedPressure}(:mixer_1, :feed, :recycle, :mix_out),
    splitter_1 = TwoSplitter{:Ratio}(:splitter_1, :c1_out, :recycle, :product),
    compressor_1 = Compressor{:PressureRatio}(:compressor_1, :mix_out, :c1_out, :electricity_1),
)

tear_points = (
    tear_point_1 = TearPoint(
        :splitter_1, :mix_1
    ),
)

function generate_graph_and_sort_topo(unordered_topo, tear_points)
    g = MetaGraph(
        DiGraph(), 
        label_type = Symbol, 
        vertex_data_type = Union{AbstractComponent, AbstractStream},
        edge_data_type = Symbol,
    )

    source_dict = Dict{Symbol, Symbol}()
    dest_dict = Dict{Symbol, Symbol}()

    for component_name in keys(unordered_topo)
        component = unordered_topo[component_name]
        
        inlet_material_streams = get_inlet_material_streams(component)
        inlet_energy_streams = get_inlet_energy_streams(component)
        outlet_material_streams = get_outlet_material_streams(component)
        outlet_energy_streams = get_outlet_energy_streams(component)

        if !isempty(outlet_material_streams)
            for stream_name in outlet_material_streams
                source_dict[stream_name] = component_name
            end
        end
        if !isempty(outlet_energy_streams)
            for stream_name in outlet_energy_streams
                source_dict[stream_name] = component_name
            end
        end

        if !isempty(inlet_material_streams)
            for stream_name in inlet_material_streams
                dest_dict[stream_name] = component_name
            end
        end
        if !isempty(inlet_energy_streams)
            for stream_name in inlet_energy_streams
                dest_dict[stream_name] = component_name
            end
        end
    end

    for component_name in keys(unordered_topo)
        component = unordered_topo[component_name]
        add_vertex!(g, component_name, component)
    end

    for component_name in keys(unordered_topo)
        component = unordered_topo[component_name]

        inlet_material_streams = get_inlet_material_streams(component)
        inlet_energy_streams = get_inlet_energy_streams(component)
        outlet_material_streams = get_outlet_material_streams(component)
        outlet_energy_streams = get_outlet_energy_streams(component)

        for stream_name in inlet_material_streams
            if haskey(source_dict, stream_name)
                neighboring_component = source_dict[stream_name]
                add_edge!(g, neighboring_component, component_name, stream_name) 
            else
                add_vertex!(g, stream_name, MaterialStream(stream_name))
                add_edge!(g, stream_name, component_name, stream_name)
            end
        end

        for stream_name in inlet_energy_streams
            if haskey(source_dict, stream_name)
                neighboring_component = source_dict[stream_name]
                add_edge!(g, neighboring_component, component_name, stream_name)
            else
                add_vertex!(g, stream_name, EnergyStream(stream_name))
                add_edge!(g, stream_name, component_name, stream_name)
            end
        end

        for stream_name in outlet_material_streams
            if !haskey(dest_dict, stream_name)
                add_vertex!(g, stream_name, MaterialStream(stream_name))
                add_edge!(g, component_name, stream_name, stream_name)
            end
        end

        for stream_name in outlet_energy_streams
            if !haskey(dest_dict, stream_name)
                add_vertex!(g, stream_name, EnergyStream(stream_name))
                add_edge!(g, component_name, stream_name, stream_name)
            end
        end
    end

    g_with_cycles = deepcopy(g)

    cycles = [[label_for(g, code) for code in cycle] for cycle in simplecycles(g_with_cycles)]
    filtered_cycles = Vector{Symbol}[]
    for cycle in cycles
        filtered_cycle = filter(x -> g[x] isa AbstractComponent, cycle)
        if !isempty(filtered_cycle)
            push!(filtered_cycles, filtered_cycle)
        end
    end

    #clean out all cycles from the main graph so topological_sort_by_dfs(g) can actually work
    for tear_point in tear_points
        edge_label = g[tear_point.component_1, tear_point.component_2]
        MetaGraphsNext.delete!(g, tear_point.component_1, tear_point.component_2) 
        #delete! must be used instead of rem_edge!()
        
        torn_out = Symbol(edge_label, "_torn_out")
        torn_in  = Symbol(edge_label, "_torn_in")
        
        # Add the product boundary to component 1
        add_vertex!(g, torn_out, MaterialStream(edge_label))
        add_edge!(g, tear_point.component_1, torn_out, edge_label)
        
        # Add the feed boundary to component 2
        add_vertex!(g, torn_in, MaterialStream(edge_label))
        add_edge!(g, torn_in, tear_point.component_2, edge_label)
    end

    topological_labels = [label_for(g, code) for code in topological_sort_by_dfs(g)]
    filter!(x -> g[x] isa AbstractComponent, topological_labels)

    ordered_topo = (; (k => unordered_topo[k] for k in topological_labels)...)

    return ordered_topo, g, g_with_cycles, filtered_cycles
end

ordered_topo, g, g_with_cycles, cycles = generate_graph_and_sort_topo(unordered_topo, tear_points)

ordered_topo

g

g_with_cycles

cycles

vertex_labels = String.(labels(g_with_cycles) |> collect)
graph_edge_labels = [string(g_with_cycles[label_for(g_with_cycles, src(e)), label_for(g_with_cycles, dst(e))]) for e in edges(g_with_cycles.graph)]
fig = Figure(size = (600, 600))
ax = Axis(fig[1, 1], title = "MetaGraphsNext + GraphMakie Example")

# Use Stress() or Spring() for layout
# Note: Graphplot uses the graphplot(graph) recipe from Graphs.jl
p = graphplot!(
    ax, 
    g_with_cycles.graph; # Extract the core SimpleGraph / DiGraph structure
    layout = Stress(),
    nlabels = string.(vertex_labels),
    nlabels_color = :black,
    node_color = :blue,
    node_size = 25,
    elabels = graph_edge_labels,
    elabels_color = :gray,
    arrow_size = 12
)

hidedecorations!(ax)
hidespines!(ax)
display(fig)