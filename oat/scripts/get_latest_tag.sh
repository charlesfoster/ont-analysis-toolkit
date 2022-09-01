#!/bin/bash
# Author: Charles Foster
# Modified from https://gist.github.com/itaysk/c023de03fe74dd3d5db336b7f9699b6b
# Purpose: get latest tag on either dockerhub or quay.io for a container

usage() {
  printf "\nUsage: bash get_latest_tag.sh --source <dockerhub / quay> --repo <repo> --image <image> [--no-alphabet]  \n"
  printf "\nReturns the latest tag for a given container image\n"
  printf "* Requires 'curl' and 'jq'\n"
  printf "* --no-alphabet = removes alphabet characters from versions before sorting, USE WITH CAUTION ONLY FOR SPECIFIC IMAGES\n\n"
}

if [[ $# -eq 0 ]] ; then
	usage
    exit 1
fi

while [ "$1" != "" ]; do
    FASTA=$1
    case $1 in
      -s | --source )         shift
                              SOURCE="$1"
                              ;;
      -r | --repo )           shift
                              REPO="$1"
                              ;;
      -i | --image )          shift
                              IMAGE="$1"
                              ;;
      -n | --no-alphabet )    shift
                              NO_ALPHABET="TRUE"
                              ;;
      -h | --help )           usage
                              exit
                              ;;
    esac
    shift
done

if [ -z "$SOURCE" ] || [ -z "$REPO" ] || [ -z "$IMAGE" ] ;
then
  echo ""
  echo "Need to provide source, repo and image."
  exit 1
fi

if [ -z "$NO_ALPHABET" ]; then
  NO_ALPHABET=FALSE
fi

if [[ ${SOURCE} == "dockerhub" ]]; then
  LATEST=$(curl -s -L --fail "https://hub.docker.com/v2/repositories/${REPO}/${IMAGE}/tags/?page_size=1000" | \
    jq .results[].name -r)
  
  if [[ $NO_ALPHABET == "TRUE" ]]; then
    FINAL=$( echo $LATEST | tr " " "\n" | egrep -v "[A-Za-z]" | sort --version-sort | tail -n 1)
  else
    FINAL=$( echo $LATEST | tr " " "\n" | egrep -v "^latest$" | sort --version-sort | tail -n 1)
  fi

elif [[ ${SOURCE} == "quay" ]]; then
  LATEST=$(curl -s -H "Authorization: Bearer XYZ" -X GET "https://quay.io/api/v1/repository/${REPO}/${IMAGE}/tag/" | \
  jq .tags[].name -r | head -1 )
  FINAL=$LATEST
  
else
  echo ""
  echo "Repo must be dockerhub or quay"
  exit 1
fi

printf $FINAL
exit 0