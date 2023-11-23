#!/bin/bash

# build rust code
# `bash build.sh`: compile from scratch.
#

set -eux

cd ..

set +ux
if [[ $(uname) == Darwin ]]; then
    echo "Apple"
    :
elif [[ $HOSTNAME = *.personal-genome.org ]]; then
    OS_VER=$(lsb_release -r | cut -f2,2 | cut -d'.' -f1,1)

    if [[ ${OS_VER} == "8" ]]; then
        :
        #echo "Use Centos7 since when compiled in Centos8, the program cannot run on Centos7."
    elif [[ ${OS_VER} == "7" ]]; then
        echo "Use Centos8. Cannot build in 7."
        exit 1
    fi

    source /bio/lmod/lmod/init/bash
    module use /nfs/data06/ricky/app/.modulefiles
    module purge

    if [[ ${OS_VER} == "8" ]]; then
        module load clang/16.0.2
    elif [[ ${OS_VER} == "7" ]]; then
        module load gcc/9.2.0
        module load clang/16.0.2
        # tmp: mv to install-script
        #https://users.rust-lang.org/t/unable-to-find-libclang-the-libclang-issue/70746
        #export LIBCLANG_PATH="/nfs/data06/ricky/app/clang7/16.0.2/build/bin"
    fi

elif [[ $HOSTNAME = *.obcx ]]; then
    source /work/00/gg57/share/yoshi-tools/lmod/lmod/init/bash
    module use /work/gg57/j29002/app/.modulefiles
    module purge
    module list
    module load gcc/7.5.0
    module load cmake/3.14.5
    module load clang/16.0.2
else
    echo "unknown host: $HOSTNAME"
    exit 1
fi
set -ux

export RUST_BACKTRACE=full

#if [ "$1" != "light" ]; then
if [[ $* != *light* ]]; then
    # to delete old unused binary
    cargo clean --release \
        --manifest-path ./projects_rust/Cargo.toml || exit 1

    cargo clean \
        --manifest-path ./projects_rust/Cargo.toml || exit 1
fi

#if [ "$1" == "no-manifest" ] || [ "$2" == "no-manifest" ]; then
if [[ $* == *no-feature* ]]; then
    export RUSTFLAGS='-C target-cpu=native'
    cargo build --release \
        --manifest-path ./projects_rust/Cargo.toml \
        -p genetics --bins \
        --no-default-features
    cargo build --release \
        --manifest-path ./projects_rust/Cargo.toml \
        -p boosting --bins \
        --no-default-features
else

    # this might speed up popcount and others
    export RUSTFLAGS='-C target-cpu=native'
    cargo build --release \
        --manifest-path ./projects_rust/Cargo.toml \
        -p genetics --bins
    cargo build --release \
        --manifest-path ./projects_rust/Cargo.toml \
        -p boosting --bins

    if [[ $* == *test* ]]; then
        cargo test \
            --manifest-path ./projects_rust/Cargo.toml \
            -p genetics
        cargo test \
            --manifest-path ./projects_rust/Cargo.toml \
            -p boosting
    fi

fi
