#!/bin/bash
remote_path=$(cat Get/remote_path)

rsync -ave ssh $remote_path/Figure/Landscapes/* $PWD/Figure/Landscapes/
