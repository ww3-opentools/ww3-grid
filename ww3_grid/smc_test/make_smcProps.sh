#!/bin/bash

ifort -r8 -mcmodel=medium smcProps.f90
mv a.out smcProps.exe
