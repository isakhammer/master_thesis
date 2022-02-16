
function run_netgen_container()  {
  image_name="netgen:latest"
  xhost +local:root
  XSOCK=/tmp/.X11-unix
  docker run -it --rm \
     -e DISPLAY=$DISPLAY \
      -v $(pwd)/:$HOME/project_thesis \
      -v $XSOCK:$XSOCK \
      -v $HOME/.ssh:/root/.ssh \
       -v $HOME/.Xauthority:/root/.Xauthority \
       --name netgen \
        --privileged \
        $image_name "$@"
}

run_netgen_container
