#!/bin/bash
set -e

dir="$(dirname $(readlink -f ${BASH_SOURCE[0]}))"
binary="${2:-$dir/renderer}"
framerate="${1:-30}"

if [[ ! -x "$binary" ]]; then
    echo "Error: Binary either doesn't exist or isn't executable."
    exit 1
fi

"$binary" | mpv -\
    --demuxer=rawvideo --demuxer-rawvideo-w=720\
    --demuxer-rawvideo-h=720 --demuxer-rawvideo-format=RGBA\
    --demuxer-rawvideo-fps="$framerate"

