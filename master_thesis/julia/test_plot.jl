using Plots

function u_ex(x)
    (x[1]^2 + x[2]^2 - 1)^2 * sin(2π * x[1]) * cos(2π * x[2])
end

function is_inside_circle(x, y, r)
    x^2 + y^2 <= r^2
end

n = 1000
x_range = range(-1, 1, length=n)
y_range = range(-1, 1, length=n)
radius = 1

u_ex_values = [is_inside_circle(x, y, radius) ? u_ex([x, y]) : NaN for x in x_range, y in y_range]

heatmap(x_range, y_range, u_ex_values, xlabel="x[1]", ylabel="x[2]", color=:viridis, colorbar=true)
