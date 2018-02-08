#!/bin/bash

WORK=/user
scripts=$WORK/'scripts'
if [ ! -d ${scripts} ]; then
        mkdir ${scripts}
else
        echo "${scripts} exists"
fi

cp -r * $scripts

