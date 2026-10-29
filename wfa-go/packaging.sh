#!/usr/bin/env sh

set -eu

# export GOEXPERIMENT=greenteagc # for go1.25
CGO_ENABLED=0 gox -os="windows darwin linux" -arch="amd64 arm64" -tags netgo -ldflags '-w -s' -asmflags '-trimpath' \
    -output "wfa-go_{{.OS}}_{{.Arch}}"

dir=binaries
rm -rf "$dir";
mkdir -p "$dir";

for f in wfa-go_*; do
    mkdir -p $dir/$f;
    mv $f $dir/$f;
    cd $dir/$f;
    mv $f $(echo $f | perl -pe 's/_[^\.]+//g');
    tar -zcf $f.tar.gz wfa-go*;
    mv *.tar.gz ../;
    cd ..;
    rm -rf $f;
    cd ..;
done;
