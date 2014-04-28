#!/bin/bash


fit1=$1
fit2=$2

total_ind=$3
ger_alvo=$4

a=0
b=0
linhas=0

let a=$total_ind*$ger_alvo
let b=2*$ger_alvo
let linhas=$a+$b

head -n $linhas $fit1 | tail -n $total_ind | sed 's/\ /\ \ \ \ \ \ \ \ \ /g' | cut -c 10-100 | sed 's/\ //g' > temporario_fit1
head -n $linhas $fit2 | tail -n $total_ind | sed 's/\ /\ \ \ \ \ \ \ \ \ /g' | cut -c 10-100 | sed 's/\ //g' > temporario_fit2

if [ -e front_"$ger_alvo".front ]; then rm front_"$ger_alvo".front; fi
touch front_"$ger_alvo".front

echo "@    title \"Fronteiras\"" >> front_"$ger_alvo".front
echo "@    subtitle \"Geracao "$ger_alvo"\"" >> front_"$ger_alvo".front
echo "@    xaxis  label \""$fit1"\"" >> front_"$ger_alvo".front
echo "@    yaxis  label \""$fit2"\"" >> front_"$ger_alvo".front
echo "@    s0 symbol 1" >> front_"$ger_alvo".front
echo "@    s0 symbol size 0.300000"  >> front_"$ger_alvo".front
echo "@    s0 symbol linewidth 3.0" >> front_"$ger_alvo".front
echo "@    s0 line type 0" >> front_"$ger_alvo".front

paste temporario_fit1 temporario_fit2 >> front_"$ger_alvo".front

rm temporario_*
