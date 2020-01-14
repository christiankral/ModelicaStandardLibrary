within Modelica.Magnetic.QuasiStatic.FluxTubes.Examples;
package BasicExamples "Educational examples"
extends Modelica.Icons.ExamplesPackage;

  model QuadraticCoreAirgap "Educational example: iron core with airgap"
    extends Modelica.Icons.Example;
    parameter Modelica.SIunits.Length l=0.1 "Outer length of iron core";
    parameter Modelica.SIunits.Length a=0.01 "Side length of square cross section";
    parameter Real mu_r=1000 "Relative permeability of core";
    parameter Modelica.SIunits.Length delta=0.001 "Length of airgap";
  parameter Integer N=500 "Number of turns of exciting coil";
    parameter Modelica.SIunits.Current I=1.5 "Maximum exciting current";
    Modelica.Magnetic.QuasiStatic.FluxTubes.Basic.ElectroMagneticConverter excitingCoil(N=N) annotation (Placement(transformation(extent={{-50,-10},{-30,10}})));
    Modelica.Magnetic.QuasiStatic.FluxTubes.Shapes.FixedShape.Cuboid leftLeg(
      nonLinearPermeability=false,
      mu_rConst=mu_r,
      l=l - a,
      a=a,
      b=a) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-30,30})));
    Modelica.Magnetic.QuasiStatic.FluxTubes.Shapes.FixedShape.Cuboid upperYoke(
      nonLinearPermeability=false,
      mu_rConst=mu_r,
      l=l - a,
      a=a,
      b=a) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=0,
          origin={0,50})));
    Modelica.Magnetic.QuasiStatic.FluxTubes.Shapes.FixedShape.Cuboid rightLeg(
      nonLinearPermeability=false,
      mu_rConst=mu_r,
      l=l - a - delta,
      a=a,
      b=a) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={30,30})));
    Modelica.Magnetic.QuasiStatic.FluxTubes.Shapes.FixedShape.Cuboid airGap(
      nonLinearPermeability=false,
      mu_rConst=1,
      l=delta,
      a=a,
      b=a) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={30,-30})));
    Modelica.Magnetic.QuasiStatic.FluxTubes.Basic.ElectroMagneticConverter measuringCoil(N=1) annotation (Placement(transformation(extent={{50,-10},{30,10}})));
    Modelica.Magnetic.QuasiStatic.FluxTubes.Shapes.FixedShape.Cuboid lowerYoke(
      nonLinearPermeability=false,
      mu_rConst=mu_r,
      l=l - a,
      a=a,
      b=a) annotation (Placement(transformation(
          extent={{10,-10},{-10,10}},
          rotation=0,
          origin={0,-50})));
    Modelica.Magnetic.QuasiStatic.FluxTubes.Basic.Ground magneticGround annotation (Placement(transformation(extent={{-40,-70},{-20,-50}})));
    Modelica.Electrical.QuasiStatic.SinglePhase.Basic.Ground electricGround1 annotation (Placement(transformation(extent={{-70,-30},{-50,-10}})));
    Modelica.Electrical.QuasiStatic.SinglePhase.Sources.VariableCurrentSource rampCurrent(
      I=I,
      duration=0.015,
      startTime=0.01) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-70,0})));
    Modelica.Magnetic.QuasiStatic.FluxTubes.Sensors.MagneticFluxSensor magFluxSensor annotation (Placement(transformation(
          extent={{10,-10},{-10,10}},
          rotation=270,
          origin={-30,-30})));
    Modelica.Electrical.QuasiStatic.SinglePhase.Basic.Ground electricGround2 annotation (Placement(transformation(extent={{50,-30},{70,-10}})));
    Modelica.Electrical.QuasiStatic.SinglePhase.Sensors.VoltageSensor voltageSensor annotation (Placement(transformation(
          extent={{-10,10},{10,-10}},
          rotation=270,
          origin={70,0})));
    Modelica.Magnetic.QuasiStatic.FluxTubes.Basic.LeakageWithCoefficient leakage(c_usefulFlux=0.9) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={0,0})));
    Modelica.Blocks.Sources.RealExpression usefulReluctance(y=1/airGap.G_m) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={20,-20})));
  equation
    connect(excitingCoil.port_p, leftLeg.port_p)
      annotation (Line(points={{-30,10},{-30,20}}, color={255,127,0}));
    connect(leftLeg.port_n, upperYoke.port_p)
      annotation (Line(points={{-30,40},{-30,50},{-10,50}}, color={255,127,0}));
    connect(upperYoke.port_n, rightLeg.port_p)
      annotation (Line(points={{10,50},{30,50},{30,40}}, color={255,127,0}));
    connect(rightLeg.port_n, measuringCoil.port_p)
      annotation (Line(points={{30,20},{30,10}}, color={255,127,0}));
    connect(measuringCoil.port_n, airGap.port_p)
      annotation (Line(points={{30,-10},{30,-20}}, color={255,127,0}));
    connect(lowerYoke.port_p, airGap.port_n)
      annotation (Line(points={{10,-50},{30,-50},{30,-40}}, color={255,127,0}));
    connect(lowerYoke.port_n, magneticGround.port)
      annotation (Line(points={{-10,-50},{-30,-50}}, color={255,127,0}));
    connect(excitingCoil.port_n, magFluxSensor.port_n)
      annotation (Line(points={{-30,-10},{-30,-20}}, color={255,127,0}));
    connect(magFluxSensor.port_p, magneticGround.port)
      annotation (Line(points={{-30,-40},{-30,-50}}, color={255,127,0}));
    connect(leakage.port_n, airGap.port_n) annotation (Line(points={{0,
            -10},{0,-10},{0,-40},{30,-40}}, color={255,127,0}));
    connect(measuringCoil.port_p, leakage.port_p)
      annotation (Line(points={{30,10},{0,10}}, color={255,127,0}));
    connect(leakage.R_mUsefulTot, usefulReluctance.y)
      annotation (Line(points={{12,0},{20,0},{20,-9}}, color={0,0,127}));
    connect(rampCurrent.pin_n, excitingCoil.pin_p) annotation (Line(points={{-70,10},{-50,10}}, color={85,170,255}));
    connect(rampCurrent.pin_p, electricGround1.pin) annotation (Line(points={{-70,-10},{-60,-10}}, color={85,170,255}));
    connect(electricGround1.pin, excitingCoil.pin_n) annotation (Line(points={{-60,-10},{-50,-10}}, color={85,170,255}));
  connect(measuringCoil.pin_p, voltageSensor.pin_p) annotation (Line(points={{50,10},{70,10}}, color={85,170,255}));
  connect(measuringCoil.pin_n, electricGround2.pin) annotation (Line(points={{50,-10},{60,-10}}, color={85,170,255}));
  connect(electricGround2.pin, voltageSensor.pin_n) annotation (Line(points={{60,-10},{70,-10}}, color={85,170,255}));
    annotation (Documentation(info="<html>
<p>
Educational example of a magnetic circuit containing an iron core and an airgap:
</p>
<p>
<img src=\"modelica://Modelica/Resources/Images/Magnetic/FluxTubes/Examples/MagneticCircuit.png\" alt=\"Magnetic circuit with iron core and airgap\">
</p>
<p>
A current ramp is applied in positive electric direction through the exciting coil, causing a rising magnetomotive force (mmf) in positive magnetic direction of the electromagnetic converter.  
The mmf in turn causes a magnetic flux through the circuit in the direction indicated by the flux sensor. 
From that magnetic flux, flux density can be calculated in every element of the magnetic circuit. Flux density is used to derive magnetic field strength.
Magnetic field strength times length of the flux line gives magnetic potential difference of each element. 
The sum of all magnetic potential differences is covered by the mmf of the exciting coil.
</p>
<p>
Using the values shown in section Parameters, the results can be validated easily by analytic calculations:
</p>
<table border=1 cellspacing=0 cellpadding=2>
<tr><th>element   </th><th>cross section</th><th>length       </th><th>rel. permeability </th><th>B                   </th><th>H                                    </th><th>mmf              </th></tr>
<tr><td>left leg  </td><td>a*a          </td><td>l - a        </td><td>&mu;<sub>r</sub>  </td><td>flux / cross section</td><td>B/(&mu;<sub>r</sub>*&mu;<sub>0</sub>)</td><td>H*length         </td></tr>
<tr><td>upper yoke</td><td>a*a          </td><td>l - a        </td><td>&mu;<sub>r</sub>  </td><td>flux / cross section</td><td>B/(&mu;<sub>r</sub>*&mu;<sub>0</sub>)</td><td>H*length         </td></tr>
<tr><td>right leg </td><td>a*a          </td><td>l - a - delta</td><td>&mu;<sub>r</sub>  </td><td>flux / cross section</td><td>B/(&mu;<sub>r</sub>*&mu;<sub>0</sub>)</td><td>H*length         </td></tr>
<tr><td>airgap    </td><td>a*a          </td><td>delta        </td><td>1                 </td><td>flux / cross section</td><td>B/&mu;<sub>0</sub>                   </td><td>H*length         </td></tr>
<tr><td>lower yoke</td><td>a*a          </td><td>l - a        </td><td>&mu;<sub>r</sub>  </td><td>flux / cross section</td><td>B/(&mu;<sub>r</sub>*&mu;<sub>0</sub>)</td><td>H*length         </td></tr>
<tr><td>total     </td><td>             </td><td>             </td><td>                  </td><td>                    </td><td>                                     </td><td>&Sigma; mmf = N*I</td></tr>
</table>
<p>
Note that since no leakage is present, the magnetic flux is the same in every element - they are connected in series. 
For calculation of the length of flux lines, a medium flux line (dashed line) is used.
</p>
<p>
Additionally, a measuring coil is placed in the airgap. 
Due to Faraday's law, the time derivative of flux causes an induced voltage both in the exciting coil (in positive direction) and in the measuring coil (in negative direction). 
Since current and therefore flux are a linear time dependent ramp, induced voltages are constant during that ramp and zero otherwise. 
Note that usage of nonlinear magnetic material would change that result due the nonlinear relationship between magnetic field strength and flux density.
</p>
<p>
Note the proper usage of electric and magnetic grounds to define zero potential.
</p>
</html>"),   experiment(StopTime=0.05, Interval=0.0001));
  end QuadraticCoreAirgap;

  model ToroidalCoreAirgap "Educational example: iron core with airgap"
    extends Modelica.Icons.Example;
    import Modelica.Constants.pi;
    parameter Modelica.SIunits.Length r=0.05 "Middle radius of iron core";
    parameter Modelica.SIunits.Length d=0.01 "Diameter of cylindrical cross section";
    parameter Modelica.SIunits.RelativePermeability mu_r=1000 "Relative permeability of core";
    parameter Modelica.SIunits.Length delta=0.001 "Length of airgap";
    parameter Modelica.SIunits.Angle alfa=(1 - delta/(2*pi*r))*2*pi "Section angle of toroidal core";
    parameter Integer N=500 "Number of exciting coil turns";
    parameter Modelica.SIunits.Current I=1.5 "Maximum exciting current";
    Modelica.Magnetic.QuasiStatic.FluxTubes.Basic.ElectroMagneticConverter excitingCoil(N=N)
      annotation (Placement(transformation(extent={{-40,-10},{-20,10}})));
    Modelica.Magnetic.QuasiStatic.FluxTubes.Shapes.FixedShape.Toroid core(
      nonLinearPermeability=false,
      mu_rConst=mu_r,
      r=r,
      d=d,
      alfa=alfa) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=0,
          origin={0,30})));
    Modelica.Magnetic.QuasiStatic.FluxTubes.Shapes.FixedShape.Toroid airGap(
      nonLinearPermeability=false,
      mu_rConst=1,
      r=r,
      d=d,
      alfa=2*pi - alfa) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=180,
          origin={0,-50})));
    Modelica.Magnetic.QuasiStatic.FluxTubes.Basic.ElectroMagneticConverter measuringCoil(N=1)
      annotation (Placement(transformation(extent={{40,-10},{20,10}})));
    Modelica.Magnetic.QuasiStatic.FluxTubes.Basic.Ground magneticGround
      annotation (Placement(transformation(extent={{-30,-70},{-10,-50}})));
    Modelica.Electrical.QuasiStatic.SinglePhase.Basic.Ground electricGround1
      annotation (Placement(transformation(extent={{-60,-30},{-40,-10}})));
    Modelica.Electrical.QuasiStatic.SinglePhase.Sources.VariableCurrentSource rampCurrent(
      I=I,
      duration=0.015,
      startTime=0.01) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-60,0})));
    Modelica.Magnetic.QuasiStatic.FluxTubes.Sensors.MagneticFluxSensor magFluxSensor
      annotation (Placement(transformation(
          extent={{10,-10},{-10,10}},
          rotation=270,
          origin={-20,-30})));
    Modelica.Electrical.QuasiStatic.SinglePhase.Basic.Ground electricGround2
      annotation (Placement(transformation(extent={{40,-30},{60,-10}})));
    Modelica.Electrical.QuasiStatic.SinglePhase.Sensors.VoltageSensor voltageSensor annotation (
        Placement(transformation(
          extent={{-10,10},{10,-10}},
          rotation=270,
          origin={60,0})));
  equation
    connect(core.port_n, measuringCoil.port_p)
      annotation (Line(points={{10,30},{20,30},{20,10}}, color={255,127,0}));
    connect(measuringCoil.port_n, airGap.port_p)
      annotation (Line(points={{20,-10},{20,-50},{10,-50}},
                                                   color={255,127,0}));
    connect(electricGround1.p, excitingCoil.n)
      annotation (Line(points={{-50,-10},{-40,-10}}, color={0,0,255}));
    connect(rampCurrent.p, electricGround1.p)
      annotation (Line(points={{-60,-10},{-50,-10}}, color={0,0,255}));
    connect(rampCurrent.n, excitingCoil.p)
      annotation (Line(points={{-60,10},{-40,10}}, color={0,0,255}));
    connect(measuringCoil.n, electricGround2.p)
      annotation (Line(points={{40,-10},{50,-10}}, color={0,0,255}));
    connect(measuringCoil.p, voltageSensor.p)
      annotation (Line(points={{40,10},{60,10}}, color={0,0,255}));
    connect(electricGround2.p, voltageSensor.n)
      annotation (Line(points={{50,-10},{60,-10}}, color={0,0,255}));
    connect(excitingCoil.port_n, magFluxSensor.port_n)
      annotation (Line(points={{-20,-10},{-20,-20}}, color={255,127,0}));
    connect(magFluxSensor.port_p, magneticGround.port)
      annotation (Line(points={{-20,-40},{-20,-50}}, color={255,127,0}));
    connect(magneticGround.port, airGap.port_n)
      annotation (Line(points={{-20,-50},{-10,-50}}, color={255,127,0}));
    connect(excitingCoil.port_p, core.port_p)
      annotation (Line(points={{-20,10},{-20,30},{-10,30}}, color={255,127,0}));
    annotation (Documentation(info="<html>
<p>
Educational example of a magnetic circuit containing a toroidal iron core with circular cross section and an airgap:
</p>
<p>
A current ramp is applied in positive electric direction through the exciting coil, causing a rising magnetomotive force (mmf) in positive magnetic direction of the electromagnetic converter.  
The mmf in turn causes a magnetic flux through the circuit in the direction indicated by the flux sensor. 
From that magnetic flux, flux density can be calculated in every element of the magnetic circuit. Flux density is used to derive magnetic field strength.
Magnetic field strength times length of the flux line gives magnetic potential difference of each element. 
The sum of all magnetic potential differences is covered by the mmf of the exciting coil.
</p>
<p>
Using the values shown in section Parameters, the results can be validated easily by analytic calculations:
</p>
<table border=1 cellspacing=0 cellpadding=2>
<tr><th>element   </th><th>cross section     </th><th>length             </th><th>rel. permeability </th><th>B                   </th><th>H                                    </th><th>mmf              </th></tr>
<tr><td>core      </td><td>d<sup>2</sup>*pi/4</td><td>r*alfa             </td><td>&mu;<sub>r</sub>  </td><td>flux / cross section</td><td>B/(&mu;<sub>r</sub>*&mu;<sub>0</sub>)</td><td>H*length         </td></tr>
<tr><td>airgap    </td><td>d<sup>2</sup>*pi/4</td><td>delta=r*(2*pi-alfa)</td><td>1</td><td>flux / cross section</td><td>B/(&mu;<sub>0</sub>)</td><td>H*delta         </td></tr>
<tr><td>total     </td><td>                  </td><td>                   </td><td>                  </td><td>                    </td><td>                                     </td><td>&Sigma; mmf = N*I</td></tr>
</table>
<p>
Note that since no leakage is present, the magnetic flux is the same in every element - they are connected in series. 
For calculation of the length of flux lines, a flux line in the middle of the toroid is used.
</p>
<p>
Additionally, a measuring coil is placed in the airgap. 
Due to Faraday's law, the time derivative of flux causes an induced voltage both in the exciting coil (in positive direction) and in the measuring coil (in negative direction). 
Since current is given as a linear-time dependent ramp, the induced voltages during that ramp are constant and otherwise zero. 
Note that usage of nonlinear magnetic material would change that result due the nonlinear relationship between magnetic field strength and flux density.
</p>
<p>
Note the proper usage of electric and magnetic grounds to define zero potential.
</p>
</html>"),   experiment(StopTime=0.05, Interval=0.0001));
  end ToroidalCoreAirgap;

  model ToroidalCoreQuadraticCrossSection
    "Educational example: iron core with airgap"
    extends Modelica.Icons.Example;
    import Modelica.Constants.pi;
    parameter Modelica.SIunits.Length r_o=0.055 "Outer radius of iron core";
    parameter Modelica.SIunits.Length r_i=0.045 "Inner radius of iron core";
    parameter Modelica.SIunits.Length l=0.01 "Length of rectangular cross section";
    parameter Modelica.SIunits.RelativePermeability mu_r=1000 "Relative permeability of core";
    parameter Modelica.SIunits.Length delta=0.001 "Length of airgap";
    parameter Modelica.SIunits.Angle alfa=(1 - delta/(2*pi*(r_o + r_i)/2))*2*pi "Section angle of toroidal core";
    parameter Integer N=500 "Number of exciting coil turns";
    parameter Modelica.SIunits.Current I=1.5 "Maximum exciting current";
    Modelica.Magnetic.QuasiStatic.FluxTubes.Basic.ElectroMagneticConverter excitingCoil(N=N)
      annotation (Placement(transformation(extent={{-40,-10},{-20,10}})));
    Modelica.Magnetic.QuasiStatic.FluxTubes.Shapes.FixedShape.HollowCylinderCircumferentialFlux core(
      nonLinearPermeability=false,
      mu_rConst=mu_r,
      l=l,
      r_i=r_i,
      r_o=r_o,
      alfa=alfa) annotation (Placement(transformation(
          extent={{-10,10},{10,-10}},
          rotation=0,
          origin={0,30})));
    Modelica.Magnetic.QuasiStatic.FluxTubes.Shapes.FixedShape.HollowCylinderCircumferentialFlux airGap(
      nonLinearPermeability=false,
      mu_rConst=1,
      l=l,
      r_i=r_i,
      r_o=r_o,
      alfa=2*pi - alfa) annotation (Placement(transformation(
          extent={{-10,10},{10,-10}},
          rotation=180,
          origin={0,-50})));
    Modelica.Magnetic.QuasiStatic.FluxTubes.Basic.ElectroMagneticConverter measuringCoil(N=1)
      annotation (Placement(transformation(extent={{40,-10},{20,10}})));
    Modelica.Magnetic.QuasiStatic.FluxTubes.Basic.Ground magneticGround
      annotation (Placement(transformation(extent={{-30,-70},{-10,-50}})));
    Modelica.Electrical.QuasiStatic.SinglePhase.Basic.Ground electricGround1
      annotation (Placement(transformation(extent={{-60,-30},{-40,-10}})));
    Modelica.Electrical.QuasiStatic.SinglePhase.Sources.VariableCurrentSource rampCurrent(
      I=I,
      duration=0.015,
      startTime=0.01) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-60,0})));
    Modelica.Magnetic.QuasiStatic.FluxTubes.Sensors.MagneticFluxSensor magFluxSensor
      annotation (Placement(transformation(
          extent={{10,-10},{-10,10}},
          rotation=270,
          origin={-20,-30})));
    Modelica.Electrical.QuasiStatic.SinglePhase.Basic.Ground electricGround2
      annotation (Placement(transformation(extent={{40,-30},{60,-10}})));
    Modelica.Electrical.QuasiStatic.SinglePhase.Sensors.VoltageSensor voltageSensor annotation (
        Placement(transformation(
          extent={{-10,10},{10,-10}},
          rotation=270,
          origin={60,0})));
  equation
    connect(core.port_n, measuringCoil.port_p)
      annotation (Line(points={{10,30},{20,30},{20,10}}, color={255,127,0}));
    connect(measuringCoil.port_n, airGap.port_p)
      annotation (Line(points={{20,-10},{20,-50},{10,-50}},
                                                   color={255,127,0}));
    connect(electricGround1.p, excitingCoil.n)
      annotation (Line(points={{-50,-10},{-40,-10}}, color={0,0,255}));
    connect(rampCurrent.p, electricGround1.p)
      annotation (Line(points={{-60,-10},{-50,-10}}, color={0,0,255}));
    connect(rampCurrent.n, excitingCoil.p)
      annotation (Line(points={{-60,10},{-40,10}}, color={0,0,255}));
    connect(measuringCoil.n, electricGround2.p)
      annotation (Line(points={{40,-10},{50,-10}}, color={0,0,255}));
    connect(measuringCoil.p, voltageSensor.p)
      annotation (Line(points={{40,10},{60,10}}, color={0,0,255}));
    connect(electricGround2.p, voltageSensor.n)
      annotation (Line(points={{50,-10},{60,-10}}, color={0,0,255}));
    connect(excitingCoil.port_n, magFluxSensor.port_n)
      annotation (Line(points={{-20,-10},{-20,-20}}, color={255,127,0}));
    connect(magFluxSensor.port_p, magneticGround.port)
      annotation (Line(points={{-20,-40},{-20,-50}}, color={255,127,0}));
    connect(magneticGround.port, airGap.port_n)
      annotation (Line(points={{-20,-50},{-10,-50}}, color={255,127,0}));
    connect(excitingCoil.port_p, core.port_p)
      annotation (Line(points={{-20,10},{-20,30},{-10,30}}, color={255,127,0}));
    annotation (Documentation(info="<html>
<p>
Educational example of a magnetic circuit containing a toroidal iron core with rectangular cross section and an airgap:
</p>
<p>
A current ramp is applied in positive electric direction through the exciting coil, causing a rising magnetomotive force (mmf) in positive magnetic direction of the electromagnetic converter.  
The mmf in turn causes a magnetic flux through the circuit in the direction indicated by the flux sensor. 
From that magnetic flux, flux density can be calculated in every element of the magnetic circuit. Flux density is used to derive magnetic field strength.
Magnetic field strength times length of the flux line gives magnetic potential difference of each element. 
The sum of all magnetic potential differences is covered by the mmf of the exciting coil.
</p>
<p>
Using the values shown in section Parameters, the results can be validated easily by analytic calculations:
</p>
<table border=1 cellspacing=0 cellpadding=2>
<tr><th>element   </th><th>cross section</th><th>length                         </th><th>rel. permeability </th><th>B                   </th><th>H                                    </th><th>mmf              </th></tr>
<tr><td>core      </td><td>(r_o - r_i)*l</td><td>(r_o + r_i)/2*alfa             </td><td>&mu;<sub>r</sub>  </td><td>flux / cross section</td><td>B/(&mu;<sub>r</sub>*&mu;<sub>0</sub>)</td><td>H*length         </td></tr>
<tr><td>airgap    </td><td>(r_o - r_i)*l</td><td>delta=(r_o + r_i)/2*(2*pi-alfa)</td><td>1</td><td>flux / cross section</td><td>B/(&mu;<sub>0</sub>)</td><td>H*delta         </td></tr>
<tr><td>total     </td><td>             </td><td>                               </td><td>                  </td><td>                    </td><td>                                     </td><td>&Sigma; mmf = N*I</td></tr>
</table>
<p>
Note that since no leakage is present, the magnetic flux is the same in every element - they are connected in series. 
For calculation of the length of flux lines, a medium flux line is used.
</p>
<p>
Additionally, a measuring coil is placed in the airgap. 
Due to Faraday's law, the time derivative of flux causes an induced voltage both in the exciting coil (in positive direction) and in the measuring coil (in negative direction). 
Since current and therefore flux are a linear time dependent ramp, induced voltages are constant during that ramp and zero otherwise. 
Note that usage of nonlinear magnetic material would change that result due the nonlinear relationship between magnetic field strength and flux density.
</p>
<p>
Note the proper usage of electric and magnetic grounds to define zero potential.
</p>
</html>"),   experiment(StopTime=0.05, Interval=0.0001));
  end ToroidalCoreQuadraticCrossSection;

annotation ();
end BasicExamples;
