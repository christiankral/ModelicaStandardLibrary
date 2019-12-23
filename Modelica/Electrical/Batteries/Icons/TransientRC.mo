within Modelica.Electrical.Batteries.Icons;
class TransientRC "Icon indictaing transient RC element usage"
  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
        Rectangle(
          origin={0.0,-25.0},
          lineColor={64,64,64},
          fillColor={255,215,136},
          fillPattern=FillPattern.Solid,
          extent={{-100.0,-75.0},{100.0,75.0}},
          radius=25.0),
        Line(
          points={{-100,0},{-10,0}},
          color={64,64,64}),
        Line(
          origin={0.0,-50.0},
          points={{-100,0},{-54,0}},
          color={64,64,64}),
        Rectangle(extent={{-54,-34},{54,-66}}, lineColor={0,0,0}),
        Line(
          origin={154,-50},
          points={{-100,0},{-54,0}},
          color={64,64,64}),
        Line(
          points={{10,0},{100,0}},
          color={64,64,64}),
        Line(points={{-10,-20},{-10,20}}, color={0,0,0}),
        Line(points={{10,-20},{10,20}}, color={0,0,0}),
        Line(points={{-54,-50},{54,-50}}, color={175,175,175}),
        Line(points={{-10,0},{10,0}}, color={175,175,175}),
        Line(points={{0,50},{0,-100}}, color={175,175,175})}),
                                Diagram(coordinateSystem(preserveAspectRatio=false)));
end TransientRC;
