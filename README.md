# qtpie-fork
Fork of QTPIE

    LDFLAGS="-lblas -llapack" ./configure --prefix=/home/andersx/opt/qtpie



How to get displaced charges: (e.g. for water.xyz)

```
andersx@nuke:~/dev/qtpie-fork/src$ ./onexyz ~/water.xyz
 Single geometry mode
 Reading in file /home/andersx/water.xyz                           
 Using QTPIE
 QTPIE Energy is  -5.4323733936084508E-002
 CHARGES -0.84723990517492864       0.42393428497415309       0.42330562020077556     
 Dipole moment (Debyes)
  0.45149768936790458        1.5674568465211423       -1.1162586068429996     
 Norm =    1.9765637876278204     
 Dipole moment (atomic units)
  0.17763282686311804       0.61668486282466417      -0.43916984851330060     
 Polarizability, analytical (atomic units)
   4.7351440089291810       0.72447552467915921       0.68233084850060022     
  0.72447552467915899        3.5184309559986087       -2.5821402083442901     
  0.68233084850059977       -2.5821402083442901        2.2163852758777476     
 Polarizability, numerical (atomic units)
   4.7351440089909991       0.72447552468442311       0.68233084850148618     
  0.72447552468442311        3.5184309560704063       -2.5821402083337635     
  0.68233084850148618       -2.5821402083337635        2.2163852759771174     
 Error due to numerical accuracy
   0.0000000000000000        0.0000000000000000        0.0000000000000000     
   0.0000000000000000        0.0000000000000000        0.0000000000000000     
   0.0000000000000000        0.0000000000000000        0.0000000000000000     
 QML          -1          -1           0 -0.84662041591714199       0.42395924753278491       0.42266116838435708     
 QML          -1           0          -1 -0.84744374772220110       0.42450678273334430       0.42293696498885680     
 QML          -1           0           0 -0.84710130558702001       0.42422271303944858       0.42287859254757137     
 QML          -1           0           1 -0.84675886345183859       0.42393864334555276       0.42282022010628584     
 QML          -1           1           0 -0.84758219525689837       0.42448617854611226       0.42309601671078612     
 QML           0          -1          -1 -0.84710145764023204       0.42395488916138507       0.42314656847884702     
 QML           0          -1           0 -0.84675901550505051       0.42367081946748908       0.42308819603756143     
 QML           0          -1           1 -0.84641657336986897       0.42338674977359331       0.42302982359627572     
 QML           0           0          -1 -0.84758234731010995       0.42421835466804880       0.42336399264206109     
 QML           0           0           0 -0.84723990517492909       0.42393428497415309       0.42330562020077606     
 QML           0           0           1 -0.84689746303974756       0.42365021528025726       0.42324724775949024     
 QML           0           1          -1 -0.84806323697998864       0.42448182017471281       0.42358141680527583     
 QML           0           1           0 -0.84772079484480756       0.42419775048081720       0.42352304436399035     
 QML           0           1           1 -0.84737835270962603       0.42391368078692121       0.42346467192270476     
 QML           1          -1           0 -0.84689761509295935       0.42338239140219375       0.42351522369076561     
 QML           1           0          -1 -0.84772094689801936       0.42392992660275353       0.42379102029526589     
 QML           1           0           0 -0.84737850476283749       0.42364585690885753       0.42373264785398002     
 QML           1           0           1 -0.84703606262765663       0.42336178721496182       0.42367427541269487     
 QML           1           1           0 -0.84785939443271618       0.42390932241552143       0.42395007201719476     
```

The lines that start with QML, has the field displacement in X Y and Z directsion, and the the list of partial charges for that displacement.

E.g, the line

   QML          -1          -1           0 -0.84662041591714199       0.42395924753278491       0.42266116838435708     

has the field displaced in E = [-de, -de, 0], with the field displacement, de = 1.0e-3 a.u.



