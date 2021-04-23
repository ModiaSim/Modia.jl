within ;
package DocuModels
  model Filter
    Modelica.Electrical.Analog.Basic.Resistor R(R=0.5)
      annotation (Placement(transformation(extent={{-50,40},{-30,60}})));
    Modelica.Electrical.Analog.Basic.Capacitor C(C=2.0) annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=-90,
          origin={-10,30})));
    Modelica.Electrical.Analog.Sources.ConstantVoltage V(V=10) annotation (
        Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=-90,
          origin={-70,30})));
  equation
    connect(V.p, R.p)
      annotation (Line(points={{-70,40},{-70,50},{-50,50}}, color={0,0,255}));
    connect(R.n, C.p)
      annotation (Line(points={{-30,50},{-10,50},{-10,40}}, color={0,0,255}));
    connect(C.n, V.n) annotation (Line(points={{-10,20},{-10,8},{-70,8},{-70,20}},
          color={0,0,255}));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)));
  end Filter;
  annotation (uses(Modelica(version="3.2.3")));
end DocuModels;
