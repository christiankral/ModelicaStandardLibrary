within Modelica.Magnetic.FluxTubes.Basic;
model ElectroMagneticConverter "Ideal electromagnetic energy conversion"

  FluxTubes.Interfaces.PositiveMagneticPort port_p "Positive magnetic port"
    annotation (Placement(transformation(extent={{90,90},{110,110}}), iconTransformation(extent={{90,90},{110,110}})));
  FluxTubes.Interfaces.NegativeMagneticPort port_n "Negative magnetic port"
    annotation (Placement(transformation(extent={{110,-110},{90,-90}}), iconTransformation(extent={{110,-110},{90,-90}})));
  Modelica.Electrical.Analog.Interfaces.PositivePin p "Positive electrical pin"
    annotation (Placement(transformation(extent={{-90,90},{-110,110}}), iconTransformation(extent={{-90,90},{-110,110}})));
  Modelica.Electrical.Analog.Interfaces.NegativePin n "Negative electrical pin"
    annotation (Placement(transformation(extent={{-110,-110},{-90,-90}}), iconTransformation(extent={{-110,-110},{-90,-90}})));
  SI.Voltage v "Voltage";
  SI.Current i(start=0, stateSelect=StateSelect.prefer) "Current";
  SI.MagneticPotentialDifference V_m "Magnetic potential difference";
  SI.MagneticFlux Phi "Magnetic flux coupled into magnetic circuit";

  parameter Real N=1 "Number of turns";

  // For information only:
  SI.MagneticFlux Psi "Flux linkage";
  SI.Inductance L_stat "Static inductance abs(Psi/i)";

protected
  Real eps=100*Modelica.Constants.eps;
equation
  v = p.v - n.v;
  0 = p.i + n.i;
  i = p.i;

  V_m = port_p.V_m - port_n.V_m;
  0 = port_p.Phi + port_n.Phi;
  Phi = port_p.Phi;

  // Converter equations
  V_m = i*N "Ampere's law";
  N*der(Phi) = -v "Faraday's law";

  // For information only
  Psi = N*Phi;
  // Use of abs() for positive results; due to Modelica sign conventions for flow into connectors
  L_stat = noEvent(if abs(i) > eps then abs(Psi/i) else abs(Psi/eps));

  annotation (
    defaultComponentName="converter",
    Diagram(coordinateSystem(
        preserveAspectRatio=false,
        extent={{-100,-100},{100,100}}), graphics={Polygon(
              points={{-134,63},{-124,60},{-134,57},{-134,63}},
              lineColor={160,160,164},
              fillColor={160,160,164},
              fillPattern=FillPattern.Solid),Line(points={{-150,60},{-125,
          60}}, color={160,160,164}),Polygon(
              points={{141,-57},{151,-60},{141,-63},{141,-57}},
              lineColor={160,160,164},
              fillColor={160,160,164},
              fillPattern=FillPattern.Solid),Line(points={{125,-60},{150,-60}},
          color={160,160,164}),Text(
              extent={{128,-56},{144,-41}},
              textColor={160,160,164},
              textString="Phi"),Text(
              extent={{128,64},{145,79}},
              textString="Phi"),Line(points={{-150,-59},{-125,-59}}, color=
          {160,160,164}),Polygon(
              points={{-140,-56},{-150,-59},{-140,-62},{-140,-56}},
              lineColor={160,160,164},
              fillColor={160,160,164},
              fillPattern=FillPattern.Solid),Text(
              extent={{-141,-56},{-124,-41}},
              textColor={160,160,164},
              textString="i"),Text(
              extent={{-150,63},{-133,78}},
              textColor={160,160,164},
              textString="i"),Line(points={{124,61},{149,61}}, color={160,160,164}),
          Polygon(
              points={{134,64},{124,61},{134,58},{134,64}},
              lineColor={160,160,164},
              fillColor={160,160,164},
              fillPattern=FillPattern.Solid)}),
    Icon(coordinateSystem(
      preserveAspectRatio=false,
      extent={{-100,-100},{100,100}}), graphics={
      Text(
        extent={{-150,150},{150,110}},
        textColor={0,0,255},
        textString="%name"),
        Line(points={{-30,100},{-90,100}},
                                         color={0,0,255}),
        Line(points={{-30,-100},{-90,-100}},
                                           color={0,0,255}),
        Line(
          points={{0,80},{-100,80}},
          color={0,0,255},
          pattern=LinePattern.Dash),
        Line(
          points={{-100,80},{-100,-80}},
          color={0,0,255},
          pattern=LinePattern.Dash),
        Line(
          points={{0,-80},{-100,-80}},
          color={0,0,255},
          pattern=LinePattern.Dash),
        Line(
          points={{100,80},{0,80}},
          color={255,127,0},
          pattern=LinePattern.Dash),
        Line(
          points={{100,-80},{0,-80}},
          color={255,127,0},
          pattern=LinePattern.Dash),
        Line(
          points={{100,80},{100,-80}},
          color={255,127,0},
          pattern=LinePattern.Dash),
        Ellipse(extent={{-4,-34},{64,34}}, lineColor={255,127,0}),
        Line(points={{30,-100},{30,-34}},color={255,127,0}),
        Line(points={{18,0},{42,0}}, color={255,127,0}),
        Line(points={{42,10},{42,-12}}, color={255,127,0}),
        Line(points={{30,34},{30,100}},color={255,127,0}),
        Line(points={{30,100},{90,100}},color={255,127,0}),
        Line(points={{30,-100},{90,-100}},
                                         color={255,127,0}),
        Text(
          extent={{-150,150},{150,110}},
          textColor={0,0,255},
          textString="%name"),
        Line(points={{18,10},{18,-12}}, color={255,127,0}),
        Line(points={{-110,30},{-110,-30}},
                                          color={0,0,255}),
        Polygon(
          points={{-110,-30},{-104,-10},{-116,-10},{-110,-30}},
          lineColor={0,0,255},
          fillColor={0,0,255},
          fillPattern=FillPattern.Solid),
        Line(points={{110,32},{110,-28}}, color={255,128,0}),
        Polygon(
          points={{110,-28},{116,-8},{104,-8},{110,-28}},
          lineColor={255,128,0},
          fillColor={255,128,0},
          fillPattern=FillPattern.Solid),
        Line(
          points={{-15,-7},{-14,-1},{-7,7},{7,7},{14,-1},{15,-7}},
          color={0,0,255},
          smooth=Smooth.Bezier,
          origin={-23,45},
          rotation=270),
        Line(
          points={{-15,-7},{-14,-1},{-7,7},{7,7},{14,-1},{15,-7}},
          color={0,0,255},
          smooth=Smooth.Bezier,
          origin={-23,15},
          rotation=270),
        Line(
          points={{-15,-7},{-14,-1},{-7,7},{7,7},{14,-1},{15,-7}},
          color={0,0,255},
          smooth=Smooth.Bezier,
          origin={-23,-15},
          rotation=270),
        Line(
          points={{-15,-7},{-14,-1},{-7,7},{7,7},{14,-1},{15,-7}},
          color={0,0,255},
          smooth=Smooth.Bezier,
          origin={-23,-45},
          rotation=270),
        Line(points={{-30,60},{-30,100}}, color={28,108,200}),
        Line(points={{-30,-100},{-30,-60}}, color={28,108,200})}),
    Documentation(info="<html>
<p>
The electromagnetic energy conversion is given by <em>Ampere</em>'s law and <em>Faraday</em>'s law respectively:
</p>

<pre>
    V_m = i * N
    N * d&Phi;/dt = -v
</pre>

<p>
V_m is the magnetomotive force that is supplied to the connected magnetic circuit, &Phi; is the magnetic flux through the associated branch of this magnetic circuit. The negative sign of the induced voltage v is due to <em>Lenz</em>'s law.
</p>

<p>
The flux linkage &Psi; and the static inductance L_stat = |&Psi;/i| are calculated for information only. Note that L_stat is set to |&Psi;/eps| if |i| &lt; eps
(= 100*Modelica.Constants.eps).
</p>
</html>"));
end ElectroMagneticConverter;
