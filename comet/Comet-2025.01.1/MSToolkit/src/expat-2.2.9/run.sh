#! /usr/bin/env bash
# Copyright (C) 2017 Expat development team
# Licensed under the MIT license

case "arm-apple-darwin24.3.0" in
*-mingw*)
    exec wine "$@"
    ;;
*)
    exec "$@"
    ;;
esac
