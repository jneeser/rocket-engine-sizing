# rocket-engine-sizing
The rocket engine sizing files is a basic librabry of rocket engine sizing tools. Tools are provided for estimating engine geometric parametrs, hot gas properties and basic liquid or gaseous injector desing. 

## Injector Sizing
Several models are provided for sizing singlke orifice injectors. Models are based on the book Liquid Rocket Thrust Chambers: Aspects of Modelling, Analysis and Design. 

## Python Wrapper for CEA
The CEA class provides a basic wrapper for [rocketcea](https://rocketcea.readthedocs.io/en/latest/quickstart.html). 

## Enging Sizing
The engine sizing tool can be used to determine basic rocket engine parametrs, such as throat area, chamber diameter and propellant mass flow. Propellants can be selected from rocketcea or from the PropLibrary.py file. Ambient conditions need to be provided. The model assumes no flow separation in the nozzle and a combustion efficiency of 1. A flow divergence efficiency is estimated.