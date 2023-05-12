namespace QSF
{
	/// @brief Base class for describing electron-nucleus interaction
	/// @param Ncharge Charge of the nucleus
	/// @param Echarge Charge of the electron
	/// @param Nsoft "Softening" parameter for the nucleus-electron interaction
	/// @param Esoft "Softening" parameter for the electron-electron interaction
	struct InteractionBase
	{
		double Ncharge;
		double Echarge;
		double Nsoft;
		double Esoft;		//,yukawaDecay;
	};

	/// @brief Basic Coulomb interaction model
	/// @note Full form of the interaction defined through operator() is valid for any number of
	/// dimensions
	/// @warning The partial interaction eliminating one or more axis out of the interaction is
	/// valid only up to 3 dimensions
	struct CoulombInteraction: InteractionBase
	{
		double NEcharge;
		static constexpr std::string_view name= "C";

		CoulombInteraction(Section& settings)
		{
			logInfo("InteractionBase init");
			inipp::get_value(settings, "Ncharge", Ncharge);
			inipp::get_value(settings, "Echarge", Echarge);
			inipp::get_value(settings, "Nsoft", Nsoft);
			inipp::get_value(settings, "Esoft", Esoft);
			NEcharge= (Ncharge * Echarge);
			logInfo("CoulombInteraction init with NEcharge=%g (%g*%g)", NEcharge, Ncharge, Echarge);
		}
		explicit CoulombInteraction(InteractionBase ib)
			: InteractionBase(ib), NEcharge(Ncharge * Echarge)
		{}

		template<typename... Args> inline double operator()(Args... coords)
		{
			// TODO: finish
			return NEcharge / sqrt(((coords * coords) + ...) + Nsoft);
		}

		inline double partial(AXIS fc, double x, double y, double z)
		{
			if(fc == AXIS::NO) return operator()(x, y, z);
			else if(fc == AXIS::Z)
				return operator()(x, y);
			else if(fc == AXIS::Y)
				return operator()(x, z);
			else if(fc == AXIS::X)
				return operator()(y, z);
			else if(fc == AXIS::YZ)
				return operator()(z);
			else if(fc == AXIS::XZ)
				return operator()(y);
			else if(fc == AXIS::XY)
				return operator()(z);
			else
				return 0.0;
		}

		inline double partial(AXIS fc, double x, double y)
		{
			if(fc == AXIS::NO) return operator()(x, y);
			else if(fc == AXIS::Y)
				return operator()(x);
			else if(fc == AXIS::X)
				return operator()(y);
			else
				return 0.0;
		}

		inline double partial(AXIS fc, double x)
		{
			if(fc == AXIS::NO) return operator()(x);
			else
				return 0.0;
		}
	};
	/// @brief Types acceptable by ReducedDimInteraction
	enum class ReducedModel
	{
		EckhardSacha,
		Eberly
	};

	/// @brief Generic reduced dimensionality model which accepts specializations ReducedModel
	/// @see http://onet.pl
	/// @tparam model Specialized
	/// @note Full form of the interaction defined through operator() is valid for any number of
	/// dimensions
	/// @warning The partial interaction eliminating one or more axis out of the interaction is
	/// valid only up to 3 dimensions
	template<ReducedModel model> struct ReducedDimInteraction: InteractionBase
	{
		double NEcharge;
		double EEcharge;
		static constexpr std::string_view name= "ES";
		explicit ReducedDimInteraction(InteractionBase p)
			: InteractionBase(p), NEcharge(Ncharge * Echarge), EEcharge(Echarge * Echarge)
		{
			logInfo(
				"ReducedDimInteraction init NEcharge=%g, EEcharge=%g, Nsoft=%g, Esoft=%g",
				NEcharge,
				EEcharge,
				Nsoft,
				Esoft);
		}

		ReducedDimInteraction(Section& settings)
		{
			logInfo("InteractionBase init");
			inipp::get_value(settings, "Ncharge", Ncharge);
			inipp::get_value(settings, "Echarge", Echarge);
			inipp::get_value(settings, "Nsoft", Nsoft);
			inipp::get_value(settings, "Esoft", Esoft);
			NEcharge= (Ncharge * Echarge);
			EEcharge= (Echarge * Echarge);
			logInfo(
				"ReducedDimInteraction init NEcharge=%g, EEcharge=%g, Nsoft=%g, Esoft=%g",
				NEcharge,
				EEcharge,
				Nsoft,
				Esoft);
		}

		inline double ee_i(double x, double y)
		{
			if constexpr(model == ReducedModel::EckhardSacha)
				return EEcharge / sqrt((x - y) * (x - y) + x * y + Esoft);
			else if constexpr(model == ReducedModel::Eberly)
				return EEcharge / sqrt((x - y) * (x - y) + Esoft);
		}
		inline double Ne_i(double x) { return NEcharge / sqrt(x * x + Nsoft); }
		inline double operator()(double x, double y, double z)
		{
			return ee_i(x, y) + ee_i(y, z) + ee_i(z, x) + Ne_i(x) + Ne_i(y) + Ne_i(z);
		}
		inline double operator()(double x, double y) { return ee_i(x, y) + Ne_i(x) + Ne_i(y); }
		inline double operator()(double x) { return Ne_i(x); }
		inline double partial(AXIS fc, double x)
		{
			if(fc == AXIS::NO) return Ne_i(x);
			else
				return 0.0;
		}
		inline double partial(AXIS fc, double x, double y)
		{
			if(fc == AXIS::NO) return operator()(x, y);
			else if(fc == AXIS::Y)
				return Ne_i(x);
			else if(fc == AXIS::X)
				return Ne_i(y);
			else
				return 0.0;
		}

		inline double partial(AXIS fc, double x, double y, double z)
		{
			if(fc == AXIS::NO) return operator()(x, y, z);
			else if(fc == AXIS::Z)
				return ee_i(x, y) + Ne_i(x) + Ne_i(y);
			else if(fc == AXIS::Y)
				return ee_i(x, z) + Ne_i(x) + Ne_i(z);
			else if(fc == AXIS::X)
				return ee_i(y, z) + Ne_i(y) + Ne_i(z);
			else if(fc == AXIS::YZ)
				return Ne_i(x);
			else if(fc == AXIS::XZ)
				return Ne_i(y);
			else if(fc == AXIS::XY)
				return Ne_i(z);
			else
				return 0.0;
		}
	};
};	 // namespace QSF