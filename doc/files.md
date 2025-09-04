# Files
  
## data

| Name                                                   | Run(s) | Data                                             | Description                                                   |
| ------------------------------------------------------ | ------ | ------------------------------------------------ | ------------------------------------------------------------- |
| g249\_all\_det\_offline\_<br>20250828\_210751.root     | –      | LosHit,<br>FrsData,<br>NeulandCal,<br>NeulanHits | Check ToF alignment                                           |
| g249\_data\_incoming\_online\_<br>20250728_183008.root | –      | FootHit,<br>LosHit                               | FOOT and LOS data after eta correction and charge calibration |

## results

| Name                                                   |  Description                                                   |
| ------------------------------------------------------ |  ------------------------------------------------------------- |
| cal/chargeId_perDet.root     | Eta vs Eenergy plot and charge correlation with LOS after eta correction applied      | 
| cal/tof_alignment_losNeuland_v2.root | ToF vs paddle number with v2 of hit parameters (WRONG)|
| cal/tof_alignment_losNeuland_LOS.root | ToF vs paddle number with v2 of hit parameters without LOS trigger|
| cal/tof_alignment_losNeuland_v1.root | ToF vs paddle number with v1 of hit parameters (WRONG)|
