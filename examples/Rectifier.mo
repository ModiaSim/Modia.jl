within ;
package Rectifier
  model Rectifier1
    extends Modelica.Icons.Example;

    Modelica.Electrical.Analog.Sources.SineVoltage sineVoltage(V=5, freqHz=1.5)
      annotation (Placement(transformation(
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
    Modelica.Electrical.Analog.Basic.Capacitor capacitor(C=0.1) annotation (
        Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=-90,
          origin={52,30})));
    Modelica.Electrical.Analog.Basic.Ground ground
      annotation (Placement(transformation(extent={{-60,-20},{-40,0}})));
    Modelica.Electrical.Analog.Ideal.IdealDiode diode
      annotation (Placement(transformation(extent={{0,40},{20,60}})));
  equation
    connect(sineVoltage.p, R1.p)
      annotation (Line(points={{-50,40},{-50,50},{-32,50}}, color={0,0,255}));
    connect(R2.n, sineVoltage.n) annotation (Line(points={{30,20},{30,6},{-50,6},{
            -50,20}}, color={0,0,255}));
    connect(capacitor.n, sineVoltage.n) annotation (Line(points={{52,20},{52,6},{-50,
            6},{-50,20}}, color={0,0,255}));
    connect(sineVoltage.n, ground.p)
      annotation (Line(points={{-50,20},{-50,0}}, color={0,0,255}));
    connect(R1.n, diode.p)
      annotation (Line(points={{-12,50},{0,50}}, color={0,0,255}));
    connect(diode.n, R2.p)
      annotation (Line(points={{20,50},{30,50},{30,40}}, color={0,0,255}));
    connect(diode.n, capacitor.p)
      annotation (Line(points={{20,50},{52,50},{52,40}}, color={0,0,255}));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)));
  end Rectifier1;

  model Rectifier2
    extends Modelica.Icons.Example;

    Modelica.Electrical.Analog.Sources.SineVoltage Vsource(V=220, freqHz=50)
      annotation (Placement(transformation(
          extent={{-10,10},{10,-10}},
          rotation=270,
          origin={-80,30})));
    Modelica.Electrical.Analog.Basic.Resistor R1(R=20)
      annotation (Placement(transformation(extent={{-62,30},{-42,50}})));
    Modelica.Electrical.Analog.Basic.Resistor R2(R=500)
                                                       annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={30,42})));
    Modelica.Electrical.Analog.Basic.Capacitor C(C=1e-4) annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=-90,
          origin={52,40})));
    Modelica.Electrical.Analog.Basic.Ground ground
      annotation (Placement(transformation(extent={{-90,-8},{-70,12}})));
    Modelica.Electrical.Analog.Ideal.IdealDiode D1 annotation (Placement(
          transformation(
          extent={{10,10},{-10,-10}},
          rotation=-90,
          origin={-10,58})));
    Modelica.Electrical.Analog.Ideal.IdealDiode D3 annotation (Placement(
          transformation(
          extent={{10,10},{-10,-10}},
          rotation=-90,
          origin={-10,0})));
    Modelica.Electrical.Analog.Ideal.IdealDiode D2 annotation (Placement(
          transformation(
          extent={{10,10},{-10,-10}},
          rotation=-90,
          origin={10,58})));
    Modelica.Electrical.Analog.Ideal.IdealDiode D4 annotation (Placement(
          transformation(
          extent={{10,10},{-10,-10}},
          rotation=-90,
          origin={10,0})));
  equation
    connect(Vsource.p, R1.p)
      annotation (Line(points={{-80,40},{-62,40}}, color={0,0,255}));
    connect(Vsource.n, ground.p)
      annotation (Line(points={{-80,20},{-80,12}}, color={0,0,255}));
    connect(R1.n, D1.p)
      annotation (Line(points={{-42,40},{-10,40},{-10,48}}, color={0,0,255}));
    connect(D3.n, D1.p)
      annotation (Line(points={{-10,10},{-10,48}}, color={0,0,255}));
    connect(D1.n, D2.n)
      annotation (Line(points={{-10,68},{10,68}}, color={0,0,255}));
    connect(D4.n, D2.p)
      annotation (Line(points={{10,10},{10,48}}, color={0,0,255}));
    connect(Vsource.n, D2.p)
      annotation (Line(points={{-80,20},{10,20},{10,48}}, color={0,0,255}));
    connect(D3.p, D4.p)
      annotation (Line(points={{-10,-10},{10,-10}}, color={0,0,255}));
    connect(D2.n, R2.p) annotation (Line(points={{10,68},{40,68},{40,56},{30,56},
            {30,52}}, color={0,0,255}));
    connect(R2.p, C.p) annotation (Line(points={{30,52},{30,56},{52,56},{52,50}},
          color={0,0,255}));
    connect(D4.p, C.n) annotation (Line(points={{10,-10},{40,-10},{40,20},{52,
            20},{52,30}}, color={0,0,255}));
    connect(R2.n, C.n) annotation (Line(points={{30,32},{30,20},{52,20},{52,30}},
          color={0,0,255}));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)));
  end Rectifier2;
  annotation (uses(Modelica(version="3.2.3")));
end Rectifier;
