@echo off
for %%e in (*.log) do del %%e
for %%f in (*.aux) do del %%f
for %%g in (*.out) do del %%g
for %%h in (*.dvi) do del %%h
for %%j in (*.blg) do del %%j
for %%j in (*.spl) do del %%j
for %%j in (*.gz) do del %%j
for %%j in (*.snm) do del %%j
for %%j in (*.toc) do del %%j
for %%j in (*.nav) do del %%j
for %%j in (*.listing) do del %%j
REM move *.pdf build/
cls