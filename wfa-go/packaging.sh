#!/usr/bin/env sh

CGO_ENABLED=0 gox -os="windows darwin linux" -arch="amd64 arm64" -tags netgo -ldflags '-w -s' -asmflags '-trimpath' \
    -output "wfa-go_{{.OS}}_{{.Arch}}"

dir=binaries
mkdir -p $dir;
rm -rf $dir/$f;

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
