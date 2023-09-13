account="stitam"
container="prophyl"
version="0.11"

# build docker image
sudo docker build -f Dockerfile -t $account/$container:$version .

# convert docker to singularity for local use
sudo singularity build $account-$container-$version.img docker-daemon://$account/$container:$version

# push to Docker Hub
# sudo docker push $account/$container:$version

