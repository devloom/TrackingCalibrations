void PHGarfield::GetElectricFieldVcm(
    double x_cm, double y_cm, double z_cm,
    double& ex_vcm, double& ey_vcm, double& ez_vcm) const
{
  const double r_cm = std::hypot(x_cm, y_cm);
  const double abs_z_cm = std::abs(z_cm);

  ex_vcm = 0.0;
  ey_vcm = 0.0;
  ez_vcm = z_cm > 0.0 ? -m_CMVoltageDefault
                      :  m_CMVoltageDefault;

  // Diagnostic static radial field.
  // k has units V/cm^2, making Er = k*r have units V/cm.
  // constexpr double kNorth = -0.05;
  constexpr double kNorth = 0.0;
  constexpr double kSouth = 0.0;

  const double k = z_cm > 0.0 ? kNorth : kSouth;
  const double test_er_vcm = k * r_cm;

  if (r_cm > 1.e-9)
  {
    ex_vcm += test_er_vcm * x_cm / r_cm;
    ey_vcm += test_er_vcm * y_cm / r_cm;
  }

  const double spaceChargeScale =
      z_cm < 0.0 ? m_spaceChargeScale_side0
                 : m_spaceChargeScale_side1;

  if (!m_erCorrection || !m_ezCorrection ||
      spaceChargeScale == 0.0)
  {
    return;
  }

  // Existing space-charge contribution continues below.
  const double delta_er_vcm =
      spaceChargeScale *
      InterpolateCorrectionVcm(m_erCorrection, r_cm, abs_z_cm);

  const double delta_ez_local_vcm =
      spaceChargeScale *
      InterpolateCorrectionVcm(m_ezCorrection, r_cm, abs_z_cm);

  if (r_cm > 0.0)
  {
    ex_vcm += delta_er_vcm * x_cm / r_cm;
    ey_vcm += delta_er_vcm * y_cm / r_cm;
  }

  ez_vcm += z_cm >= 0.0
                ? delta_ez_local_vcm
                : -delta_ez_local_vcm;
}