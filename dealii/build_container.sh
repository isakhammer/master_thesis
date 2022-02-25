
function build_container()  {
  image_name="dealii_project_thesis:latest"
  docker build . -t $image_name
}

build_container
