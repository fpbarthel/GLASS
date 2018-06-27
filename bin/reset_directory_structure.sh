#!/bin/bash
rm benchmarks/*/*
rmdir benchmarks/*
rmdir benchmarks

rm logs/*/*
rmdir logs/*
rmdir logs

mkdir -p logs/drmaa

rm results/*/*/*/*
rm results/*/*/*
rm results/*/*
rm results/*

rmdir results/*/*/*
rmdir results/*/*
rmdir results/*
rmdir results