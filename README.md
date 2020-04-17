# RBS-callibration
Rutherford backscattering spectrometry and related techniques have long been used to determine the element depth profiles in films a few nm to a few microns thick. An ion beam with typical energy of several MeV is sent to the sample, penetrates through the sample, and loses its kinetic energy “continuously” (the so-called Continuous Slowing-Down Approximation, CSDA) in a well-known manner. Along the ion trajectories, there is a chance for collisions with target nuclei. In RBS, back-scattered ions of mass M1 are detected. The scattering bodies are the target nuclei of mass M2. Assuming the elastic scattering and applying the momentum-conservation and the energy-conservation laws, the energy of the back-scattered ions, E, is unambiguously determined by the mass of a scattering target nucleus, M2 (for a defined scattering angle, θ):
	E = E0*K	(4.1)

where
	K=  (M_1^2)/(M_1+ M_2 )^2  {〖cos⁡〖θ ± 〗 [(M_2/M_1 )^2- sin^2 θ]〗^(1/2) }^2	(4.2)
 
Here K, the kinematic factor, characterises the target atom, E0 is the energy of the incident ion immediately before the collision and 𝜃 is the scattering angle. RBS (Rutherford backscattering spectrometry) detector converts the energy of the detected particles to a charge pulse and a charge-sensitive pre-amplifier generates a voltage pulse for each incident particle. The result of an RBS measurement is the yield (number of pulses, counts) in individual channels. The RBS spectra carry information about the elemental composition of the sample including the concentration depth-distributions. The accelerated ions are in most cases 4He+. A solid state detector converts the energy of the detected particles to a charge pulse and a pre-amplifier generates a voltage pulse for each detected particle. The height of the pulse, which is proportional to the energy of the detected particle, is converted to a digital output and displayed by a multi-channel analyser (MCA). [1,2] There is a linear dependence between the energy and the channel-number, ch:
	E = k*ch+q	(4.3)

where k [keV/channel] and q [keV] are the energy calibration coefficients. The aim of this project was software development for exact calculation of the energy per channel, k, and detector offset, q. These parameters are necessary for RBS spectrum analysis.
Calibration procedure
An in-house prepared sample with dimensions of 1x1x0.1 cm was used as a reference. On a C substrate, patches of Si, Ni, and Au are deposited whereby the substrate is not completely covered. The dimensions of the Si, Ni, and Au patches are 40x40 µm, 25x25 µm, and 12x12 µm. The reference was irradiated with 1700 keV 4He+ ions with the total charge of 20 µC and the RBS spectrum was measured at scattering angle 𝜃 =170°. The RBS spectrum is shown in figure 2a. The three peaks in range 300-900 channels with height of about 800 counts represent the signal from the Si (400-550), Ni (600-750), and Au (650-900) patches. The large and high peak on the left hand side represents the signal from C (0-300). The high narrow peak in the range 900-1050 is an artificially introduced charge measurement on the sample. The peak edges on the right side represent back-scattered ions from the surfaces of C and the Si, Ni, and Au patches. Their exact positions are determined by derivative of the RBS spectrum. In Figure 2c, derivative of the RBS spectrum in range 210-260 channels (C edge) is depicted by the blue curve. A spline of the derivative (green curve) is fitted by a Gaussian function (red curve), from which position of its centre is determined. In this way, software calculates the Si (Fig. 2d), Ni (Fig. 2e), and Au (Fig. 2f) edges. The energies (Ec=429.484 keV, ESi=961.744 keV, ENi=1296.28 keV, and EAu =1568.23 keV) of the particles scattered from individual element surfaces calculated according to Eq. (1) are depicted in Figure 2b, where E0 = 1700 keV, M1 = 4.0026 amu, M2C = 12.0107 amu (amu = 1.66053e-27 kg). In accordance with Eq. (2), the calibration coefficients k =1.79 keV/channel and q=17.0341 keV were determined by linear fit.