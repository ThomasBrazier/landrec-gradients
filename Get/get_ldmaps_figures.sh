#!/bin/bash
remote_path=$(cat Get/remote_path)

rm $PWD/Figure/Landscapes/Rho_sliding100kb/*
rm $PWD/Figure/Landscapes/Rho/*
rm $PWD/Figure/Landscapes/Hotspot_peak_rate/*
rm $PWD/Figure/Landscapes/Hotspot_count/*

scp $remote_path/Figure/Landscapes/Rho_sliding100kb/* $PWD/Figure/Landscapes/Rho_sliding100kb/
scp $remote_path/Figure/Landscapes/Rho/* $PWD/Figure/Landscapes/Rho/
scp $remote_path/Figure/Landscapes/Hotspot_peak_rate/* $PWD/Figure/Landscapes/Hotspot_peak_rate/
scp $remote_path/Figure/Landscapes/Hotspot_count/* $PWD/Figure/Landscapes/Hotspot_count/

