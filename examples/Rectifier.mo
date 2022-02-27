within ;
package Rectifier
  model Rectifier1
    extends Modelica.Icons.Example;

    Modelica.Electrical.Analog.Sources.SineVoltage V(V=5, f=1.5) annotation (
        Placement(transformation(
          extent={{-10,10},{10,-10}},
          rotation=270,
          origin={-50,30})));
    Modelica.Electrical.Analog.Basic.Resistor R1(R=1)
      annotation (Placement(transformation(extent={{-32,40},{-12,60}})));
    Modelica.Electrical.Analog.Basic.Resistor R2(R=10) annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={30,30})));
    Modelica.Electrical.Analog.Basic.Capacitor C(C=0.1) annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=-90,
          origin={52,30})));
    Modelica.Electrical.Analog.Basic.Ground ground
      annotation (Placement(transformation(extent={{-60,-20},{-40,0}})));
    Modelica.Electrical.Analog.Ideal.IdealDiode D
      annotation (Placement(transformation(extent={{0,40},{20,60}})));
  equation
    connect(V.p, R1.p)
      annotation (Line(points={{-50,40},{-50,50},{-32,50}}, color={0,0,255}));
    connect(R2.n, V.n) annotation (Line(points={{30,20},{30,6},{-50,6},{-50,20}},
          color={0,0,255}));
    connect(C.n, V.n) annotation (Line(points={{52,20},{52,6},{-50,6},{-50,20}},
          color={0,0,255}));
    connect(V.n, ground.p)
      annotation (Line(points={{-50,20},{-50,0}}, color={0,0,255}));
    connect(R1.n, D.p)
      annotation (Line(points={{-12,50},{0,50}}, color={0,0,255}));
    connect(D.n, R2.p)
      annotation (Line(points={{20,50},{30,50},{30,40}}, color={0,0,255}));
    connect(D.n, C.p)
      annotation (Line(points={{20,50},{52,50},{52,40}}, color={0,0,255}));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)),
      experiment(StopTime=3, __Dymola_Algorithm="Dassl"));
  end Rectifier1;

  model Rectifier2
    extends Modelica.Icons.Example;

    Modelica.Electrical.Analog.Sources.SineVoltage V(V=220, f=50) annotation (
        Placement(transformation(
          extent={{-10,10},{10,-10}},
          rotation=270,
          origin={-60,10})));
    Modelica.Electrical.Analog.Basic.Resistor R1(R=20)
      annotation (Placement(transformation(extent={{-42,10},{-22,30}})));
    Modelica.Electrical.Analog.Basic.Resistor R2(R=500)
                                                       annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={30,18})));
    Modelica.Electrical.Analog.Basic.Capacitor C(C=1e-4) annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=-90,
          origin={52,18})));
    Modelica.Electrical.Analog.Basic.Ground ground
      annotation (Placement(transformation(extent={{-70,-28},{-50,-8}})));
    Modelica.Electrical.Analog.Ideal.IdealDiode D1 annotation (Placement(
          transformation(
          extent={{10,10},{-10,-10}},
          rotation=-90,
          origin={-10,38})));
    Modelica.Electrical.Analog.Ideal.IdealDiode D3 annotation (Placement(
          transformation(
          extent={{10,10},{-10,-10}},
          rotation=-90,
          origin={-10,-20})));
    Modelica.Electrical.Analog.Ideal.IdealDiode D2 annotation (Placement(
          transformation(
          extent={{10,10},{-10,-10}},
          rotation=-90,
          origin={10,38})));
    Modelica.Electrical.Analog.Ideal.IdealDiode D4 annotation (Placement(
          transformation(
          extent={{10,10},{-10,-10}},
          rotation=-90,
          origin={10,-20})));
  equation
    connect(V.p, R1.p)
      annotation (Line(points={{-60,20},{-42,20}}, color={0,0,255}));
    connect(V.n, ground.p)
      annotation (Line(points={{-60,0},{-60,-8}}, color={0,0,255}));
    connect(R1.n, D1.p)
      annotation (Line(points={{-22,20},{-10,20},{-10,28}}, color={0,0,255}));
    connect(D3.n, D1.p)
      annotation (Line(points={{-10,-10},{-10,28}},color={0,0,255}));
    connect(D1.n, D2.n)
      annotation (Line(points={{-10,48},{10,48}}, color={0,0,255}));
    connect(D4.n, D2.p)
      annotation (Line(points={{10,-10},{10,28}},color={0,0,255}));
    connect(V.n, D2.p)
      annotation (Line(points={{-60,0},{10,0},{10,28}}, color={0,0,255}));
    connect(D3.p, D4.p)
      annotation (Line(points={{-10,-30},{10,-30}}, color={0,0,255}));
    connect(D2.n, R2.p) annotation (Line(points={{10,48},{40,48},{40,36},{30,36},
            {30,28}}, color={0,0,255}));
    connect(R2.p, C.p) annotation (Line(points={{30,28},{30,36},{52,36},{52,28}},
          color={0,0,255}));
    connect(D4.p, C.n) annotation (Line(points={{10,-30},{40,-30},{40,0},{52,0},
            {52,8}},      color={0,0,255}));
    connect(R2.n, C.n) annotation (Line(points={{30,8},{30,0},{52,0},{52,8}},
          color={0,0,255}));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)),
      experiment(StopTime=0.1, __Dymola_Algorithm="Dassl"));
  end Rectifier2;
  annotation (uses(Modelica(version="4.0.0")),
    version="1",
    conversion(from(version="", script=
            "modelica://Rectifier/ConvertFromRectifier_.mos")));
end Rectifier;
