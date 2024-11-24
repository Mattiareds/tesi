MODULE STRUCTURES

  Integer, Parameter :: nummax=2000
  Integer :: vicini(nummax,nummax)
  Integer :: num_vicini(nummax)
  Integer :: pv_matrix(nummax,nummax)
  Integer :: legami_specie1,legami_specie2,legami_misti
  Integer :: ordine_gruppo
  Real*8 :: sig555,sig322,sig433,sig321,sig331,sig422,sig421,sig333,sig554,sig666 
  Real*8 :: sig533,sig211,sig544,sig200,sig311,sig300
  Real*8  :: nn(3)
  Real*8  :: dist_matrix(nummax,nummax)
  Character*7 :: gruppo
  Character*3 :: caratt
END MODULE STRUCTURES
