
function build_container()  {
  image_name="netgen:latest"
  docker build . -t $image_name
}

build_container
