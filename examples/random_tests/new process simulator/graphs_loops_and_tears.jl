using Graphs
using MetaGraphsNext

abstract type AbstractComponent{CalcMode} end

struct Compressor{CalcMode} <: AbstractComponent{CalcMode}
    id::Symbol
    stream_in::Symbol
    stream_out::Symbol
    energy_in::Symbol
end

struct TwoMixer{CalcMode} <: AbstractComponent{CalcMode}
    id::Symbol
    stream_in_1::Symbol
    stream_in_2::Symbol
    stream_out::Symbol
end

struct TwoSplitter{CalcMode} <: AbstractComponent{CalcMode}
    id::Symbol
    stream_in::Symbol
    stream_out_1::Symbol
    stream_out_2::Symbol
end

abstract type AbstractStream end
struct MaterialStream <: AbstractStream end
struct EnergyStream <: AbstractStream end

g = MetaGraph(
    DiGraph(), 
    Label = Symbol, 
    VertexData = AbstractComponent{CalcMode}, 
    EdgeData = , 
    GraphData = Nothing
)

component_names = (
    :mixer,
    :compressor,
    :splitter
)

component_types = (
    mixer = TwoMixer{:AveragedPressure},
    compressor = Compressor{:PressureRatio},
    splitter = TwoSplitter{:Ratio}
)

material_stream_names = (
    :feed,
    :mixer_outlet,
    :compressor_outlet,
    :product,
    :recycle_1,
    :recycle_2
)

energy_stream_names = (
    :electricity_1,
)

all_names = (component_names..., material_stream_names..., energy_stream_names...)

add_vertices!(g, length(all_names))

idxs = NamedTuple{all_names}(1:length(all_names))

add_edge!(g, idxs.feed, idxs.mixer)
add_edge!(g, idxs.mixer, idxs.mixer_outlet)
add_edge!(g, idxs.mixer_outlet, idxs.compressor)
add_edge!(g, idxs.electricity_1, idxs.compressor)
add_edge!(g, idxs.compressor, idxs.compressor_outlet)
add_edge!(g, idxs.compressor_outlet, idxs.splitter)
add_edge!(g, idxs.splitter, idxs.product)
add_edge!(g, idxs.splitter, idxs.recycle_1)
add_edge!(g, idxs.recycle_1, idxs.recycle_2)
add_edge!(g, idxs.recycle_2, idxs.mixer)

function get_names_in_cycle(g)
    if is_cyclic(g)
        cycles = simplecycles(g)
        names_in_cycle = []
        for cycle in cycles
            push!(names_in_cycle, all_names[cycle])
        end
        return names_in_cycle
    end
end

torn_graph = deepcopy(g)
rem_edge!(torn_graph, idxs.recycle_2, idxs.mixer)

get_names_in_cycle(g)
topological_sort_by_dfs(torn_graph)

filter(x -> all_names[x] ∈ component_names, topological_sort_by_dfs(torn_graph))
filter(x -> x ∈ component_names, all_names[topological_sort_by_dfs(torn_graph)])

inneighbors(g, idxs.mixer)
outneighbors(g, idxs.mixer)

#somehow we have to figure out a way to 


function map_component(g, topo, idxs, ::Val{Compressor{CalcMode}}, component_name, all_names) where CalcMode
    in_neighbors = inneighbors(g, idxs[component_name])
    out_neighbors = outneighbors(g, idxs[component_name])

    material_stream_in = filter(x -> x ∈ material_stream_names, all_names[in_neighbors])[1]
    material_stream_out = filter(x -> x ∈ material_stream_names, all_names[out_neighbors])[1]
    energy_stream_in = filter(x -> x ∈ energy_stream_names, all_names[in_neighbors])[1]

    comp = Compressor{CalcMode}(component_name, material_stream_in, material_stream_out, energy_stream_in)
    return merge(topo, NamedTuple{(component_name,)}((comp,)))
end

function map_component(g, topo, idxs, ::Val{TwoMixer{CalcMode}}, component_name, all_names) where CalcMode
    in_neighbors = inneighbors(g, idxs[component_name])
    out_neighbors = outneighbors(g, idxs[component_name])

    material_stream_in_1 = filter(x -> x ∈ material_stream_names, all_names[in_neighbors])[1]
    material_stream_in_2 = filter(x -> x ∈ material_stream_names, all_names[in_neighbors])[2]
    material_stream_out = filter(x -> x ∈ material_stream_names, all_names[out_neighbors])[1]

    comp = TwoMixer{CalcMode}(component_name, material_stream_in_1, material_stream_in_2, material_stream_out)
    return merge(topo, NamedTuple{(component_name,)}((comp,)))
end

function map_component(g, topo, idxs, ::Val{TwoSplitter{CalcMode}}, component_name, all_names) where CalcMode
    in_neighbors = inneighbors(g, idxs[component_name])
    out_neighbors = outneighbors(g, idxs[component_name])

    material_stream_in = filter(x -> x ∈ material_stream_names, all_names[in_neighbors])[1]
    material_stream_out_1 = filter(x -> x ∈ material_stream_names, all_names[out_neighbors])[1]
    material_stream_out_2 = filter(x -> x ∈ material_stream_names, all_names[out_neighbors])[2]

    comp = TwoSplitter{CalcMode}(component_name, material_stream_in, material_stream_out_1, material_stream_out_2)
    return merge(topo, NamedTuple{(component_name,)}((comp,)))
end

function map_components!(g, idxs, ordered_component_names, component_types, all_names)
    topo = (;)

    for component_name in ordered_component_names
        component_type = component_types[component_name]
        topo = map_component(g, topo, idxs, Val(component_type), component_name, all_names)
    end

    return topo
end

ordered_component_names = filter(x -> x ∈ component_names, all_names[topological_sort_by_dfs(torn_graph)])

topo = map_components!(g, idxs, ordered_component_names, component_types, all_names)

#example topo
unsorted_topo = (
    mix_1 = TwoMixer{:AveragedPressure}(:mixer_1, :feed, :recycle, :mix_out),
    compressor_1 = Compressor{:PressureRatio}(:compressor_1, :mix_out, :c1_out, :electricity_1),
    splitter_1 = TwoSplitter{:Ratio}(:splitter_1, :c1_out, :recycle, :product)
)


function generate_graph_and_sort_topo(unsorted_topo)
    g = MetaGraph(
        DiGraph(), 
        Label = Symbol, 
        VertexData = AbstractComponent{CalcMode}, 
        EdgeData = AbstractStream,
    )
    
        
end