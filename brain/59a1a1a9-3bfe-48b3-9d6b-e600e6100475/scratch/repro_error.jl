using LsqFit

# Dummy data similar to what the user might have
# Timestamps > 540s, maybe up to 1000s
x = collect(550.0:10.0:1000.0)
# Temperatures around 300-400 K
y = 400.0 .- (x .- 550.0) .* 0.1 # Some linear-ish decay

model(x, p) = p[1] .* exp.(p[2] .* x)
p0 = [1.0, 1.0]

try
    fit = curve_fit(model, x, y, p0)
    println("Fit successful")
catch e
    println("Error caught: ", e)
end
