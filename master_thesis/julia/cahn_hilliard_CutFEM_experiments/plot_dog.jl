# using Pkg
# Pkg.add("Interpolations")
# Pkg.add("Images")
# Pkg.add("ImageMagick")
# Pkg.add("Plots")

using Images, Interpolations, ImageMagick, Plots

function load_jpg()
    # Load the image
    img = load("dog.jpg") # Put the path to your image file here

    # Convert to grayscale to simplify the image to a 2D array
    img_gray = Gray.(img)

    # Get the image data as an array
    array = float.(channelview(img_gray))

    # Create ranges for x and y
    xs = LinRange(-1, 1, size(array, 2)) # x corresponds to width (columns)
    ys = LinRange(-1, 1, size(array, 1)) # y corresponds to height (rows)

    # Create a grid for interpolation
    grid = (collect(ys), collect(xs))

    # Create an interpolation function
    f = extrapolate(interpolate(grid, array, Gridded(Linear())), 0)

end

f= load_jpg()
# compute the interpolated image values at each coordinate
zs = [f(y, x) for y in ys, x in xs]

# plot the heatmap
heatmap(xs, ys, zs, aspect_ratio=1)

