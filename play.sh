#!/bin/bash
set -e

dir="$(dirname $(readlink -f ${BASH_SOURCE[0]}))"
binary="$dir/${1:-tetra}"
framerate="${2:-30}"

if [[ ! -x "$binary" ]]; then
    echo "Error: Binary either doesn't exist or isn't executable."
    exit 1
fi

w="720"
h="720"
size="$((w * h * 4))"
"$binary" | mpv -\
    --demuxer=rawvideo --demuxer-rawvideo-w="$w"\
    --demuxer-rawvideo-h="$h" --demuxer-rawvideo-format=RGBA\
    --demuxer-rawvideo-size="$size" --demuxer-rawvideo-fps="$framerate"

