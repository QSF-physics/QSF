namespace QSF
{
	// Gener
	struct InteractionBase
	{
		double Ncharge;
		double Echarge;
		double Nsoft;
		double Esoft;//,yukawaDecay;
	};

	struct CoulombInteraction : InteractionBase
	{
		double NEcharge;
		static constexpr std::string_view name = "C";

		CoulombInteraction(Section& settings)
		{
			logInfo("InteractionBase init");
			inipp::get_value(settings, "Ncharge", Ncharge);
			inipp::get_value(settings, "Echarge", Echarge);
			inipp::get_value(settings, "Nsoft", Nsoft);
			inipp::get_value(settings, "Esoft", Esoft);
			NEcharge = (Ncharge * Echarge);
			logInfo("CoulombInteraction init with NEcharge=%g (%g*%g)", NEcharge, Ncharge, Echarge);
		}
		explicit CoulombInteraction(InteractionBase ib) :
			InteractionBase(ib), NEcharge(Ncharge* Echarge) {}


		template <typename ... Args>
		inline double operator()(Args...coords)
		{
			//TODO: finish
			return NEcharge / sqrt(((coords * coords) + ...) + Nsoft);
		}


		inline double partial(AXIS fc, double x, double y, double z)
		{
			if (fc == AXIS::NO)
				return operator()(x, y, z);
			else if (fc == AXIS::Z)
				return operator()(x, y);
			else if (fc == AXIS::Y)
				return operator()(x, z);
			else if (fc == AXIS::X)
				return operator()(y, z);
			else if (fc == AXIS::YZ)
				return operator()(z);
			else if (fc == AXIS::XZ)
				return operator()(y);
			else if (fc == AXIS::XY)
				return operator()(z);
			else return 0.0;
		}

		inline double partial(AXIS fc, double x, double y)
		{
			if (fc == AXIS::NO)
				return operator()(x, y);
			else if (fc == AXIS::Y)
				return operator()(x);
			else if (fc == AXIS::X)
				return operator()(y);
			else return 0.0;
		}

		inline double partial(AXIS fc, double x)
		{
			if (fc == AXIS::NO)
				return operator()(x);
			else return 0.0;
		}
	};

	struct EckhardtSachaInteraction : InteractionBase
	{
		double NEcharge;
		double EEcharge;
		static constexpr std::string_view name = "ES";
		explicit EckhardtSachaInteraction(InteractionBase p) :
			InteractionBase(p),
			NEcharge(Ncharge* Echarge),
			EEcharge(Echarge* Echarge)
		{
			logInfo("EckhardtSachaInteraction init NEcharge=%g, EEcharge=%g", NEcharge, EEcharge);
		}

		EckhardtSachaInteraction(Section& settings)
		{
			logInfo("InteractionBase init");
			inipp::get_value(settings, "Ncharge", Ncharge);
			inipp::get_value(settings, "Echarge", Echarge);
			inipp::get_value(settings, "Nsoft", Nsoft);
			inipp::get_value(settings, "Esoft", Esoft);
			NEcharge = (Ncharge * Echarge);
			EEcharge = (Echarge * Echarge);
			logInfo("EckhardtSachaInteraction init NEcharge=%g, EEcharge=%g", NEcharge, EEcharge);
		}

		inline double ee_EckhardtSacha(double x, double y)
		{
			return EEcharge / sqrt((x - y) * (x - y) + x * y + Esoft);
		}
		inline double Ne_EckhardtSacha(double x)
		{
			return NEcharge / sqrt(x * x + Nsoft);
		}
		inline double operator()(double x, double y, double z)
		{
			return ee_EckhardtSacha(x, y) + ee_EckhardtSacha(y, z) + ee_EckhardtSacha(z, x) + Ne_EckhardtSacha(x) + Ne_EckhardtSacha(y) + Ne_EckhardtSacha(z);
		}
		inline double operator()(double x, double y)
		{
			return ee_EckhardtSacha(x, y) + Ne_EckhardtSacha(x) + Ne_EckhardtSacha(y);
		}
		inline double operator()(double x)
		{
			return Ne_EckhardtSacha(x);
		}
		inline double partial(AXIS fc, double x)
		{
			if (fc == AXIS::NO)
				return Ne_EckhardtSacha(x);
			else return 0.0;
		}
		inline double partial(AXIS fc, double x, double y)
		{
			if (fc == AXIS::NO)
				return operator()(x, y);
			else if (fc == AXIS::Y)
				return Ne_EckhardtSacha(x);
			else if (fc == AXIS::X)
				return Ne_EckhardtSacha(y);
			else return 0.0;
		}

		inline double partial(AXIS fc, double x, double y, double z)
		{
			if (fc == AXIS::NO)
				return operator()(x, y, z);
			else if (fc == AXIS::Z)
				return ee_EckhardtSacha(x, y) + Ne_EckhardtSacha(x) + Ne_EckhardtSacha(y);
			else if (fc == AXIS::Y)
				return ee_EckhardtSacha(x, z) + Ne_EckhardtSacha(x) + Ne_EckhardtSacha(z);
			else if (fc == AXIS::X)
				return ee_EckhardtSacha(y, z) + Ne_EckhardtSacha(y) + Ne_EckhardtSacha(z);
			else if (fc == AXIS::YZ)
				return Ne_EckhardtSacha(x);
			else if (fc == AXIS::XZ)
				return Ne_EckhardtSacha(y);
			else if (fc == AXIS::XY)
				return Ne_EckhardtSacha(z);
			else return 0.0;
		}
	};
};