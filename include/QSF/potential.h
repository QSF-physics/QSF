struct InteractionBase
{
	/* If you think that a generic field parameter is lacking add it here */
	double Ncharge;
	double Echarge;
	double Nsoft;
	double Esoft;//,yukawaDecay;
};


// template <size_t ID, typename ... Inter>
// struct Potential
// {
// 	using type = Potential;
// 	inline static Potential instance;
// 	int count = sizeof...(Inter);

// 	InteractionBase* inters[sizeof...(Inter)]{};
// 	static constexpr std::string_view name = join_v<Inter::name...>;

// 	template <AXIS fc, typename ... Args>
// 	inline double operator()(Args...coords)
// 	{
// 		// size_t i = 0;
// 		// double val = 0;
// 		// ([&] {
// 		// 	val += static_cast<Inter*>(inters[i])->operator() < fc > (coords...) + ...);
// 		// 	i++;
// 		// }(), ...);
// 		// return val;
// 		return (static_cast<Inter*>(inters[Index_v<Inter, Inter...>])->template operator() < fc > (coords...) + ...);
// 	}
// };


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
		logInfo("CoulombInteraction init");
		NEcharge = (Ncharge * Echarge);
	}
	explicit CoulombInteraction(InteractionBase ib) :
		InteractionBase(ib), NEcharge(Ncharge* Echarge) {}


	template <typename ... Args>
	inline double operator()(Args...coords)
	{
		//TODO: finish
		return NEcharge / sqrt(((coords * coords) + ...) + Nsoft);
	}
};

struct EckhardtSachaInteraction : InteractionBase
{
	double NEcharge;
	double EEcharge;
	static constexpr std::string_view name = "ES";
	EckhardtSachaInteraction(InteractionBase p) :
		InteractionBase(p), NEcharge(Ncharge* Echarge),
		EEcharge(Echarge* Echarge) {}

	inline double ee_EckhardtSacha(double x, double y)
	{
		return EEcharge / sqrt((x - y) * (x - y) + x * y + Esoft);
	}
	inline double Ne_EckhardtSacha(double x)
	{
		return NEcharge / sqrt(x * x + Nsoft);
	}

	template <AXIS fc> inline double operator()(double x)
	{
		return Ne_EckhardtSacha(x);
	}
	template <AXIS fc> inline double operator()(double x, double y)
	{
		if (fc == AXIS::Y)
			return Ne_EckhardtSacha(x);
		if (fc == AXIS::X)
			return Ne_EckhardtSacha(y);
		else
			return ee_EckhardtSacha(x, y) + Ne_EckhardtSacha(x) + Ne_EckhardtSacha(y);
	}
	template <AXIS fc>
	inline double operator()(double x, double y, double z)
	{
		if (fc == AXIS::Z)
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
		else return ee_EckhardtSacha(x, y) + ee_EckhardtSacha(y, z) + ee_EckhardtSacha(z, x) + Ne_EckhardtSacha(x) + Ne_EckhardtSacha(y) + Ne_EckhardtSacha(z);
	}
};
