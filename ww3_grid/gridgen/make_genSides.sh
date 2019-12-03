#!/bin/bash

ifort -r8 genSides.f90
mv a.out genSides.exe
