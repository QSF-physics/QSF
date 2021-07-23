
// template <typename ... Args> struct _CO_AVG : COMPUTATION <double, Args...>
// {


	// static constexpr std::string_view name = "AVG_";
	// template<REP R, OPTIMS opt, typename Operator> inline double calc() const noexcept
	// {
		// Operator op;
		// double result = 0.0;
		// if constexpr (is_same_v<Operator, IdentityOperator<REP::X>>
		// 			  || is_same_v<Operator, IdentityOperator<REP::P>>
		// 			  || is_same_v<Operator, IdentityOperator<REP::X | REP::P>>)
		// {
		// 	for (i = 0; i < m_l; i++) result += norm(WF::psi[i]);
		// }
		// else
		// {
		// 	for (ind i = 0; i < n0_l; i++)
		// 	{
		// 		if constexpr (DIM == 1)
		// 			result += norm(WF::psi[i])
		// 			* op.template operator() < R, opt > (n0_o + i);

		// 		else
		// 		{
		// 			readInd0 = i * n;
		// 			for (ind j = 0; j < n; j++)
		// 			{
		// 				readInd1 = i * n + j;
		// 				if constexpr (DIM == 2)
		// 				{
		// 					if constexpr (R == REP::X)
		// 						result += norm(WF::psi[readInd1])
		// 						* op.template operator() < R, opt > (n0_o + i, j);
		// 					else
		// 						result += norm(WF::psi[readInd1])
		// 						* op.template operator() < R, opt > (j, n0_o + i);
		// 				}
		// 				else
		// 				{
		// 					readInd1 = readInd1 * n;
		// 					for (ind k = 0; k < n; k++)
		// 					{
		// 						if constexpr (R == REP::X)
		// 							result += norm(WF::psi[readInd1 + k])
		// 							* op.template operator() < R, opt > (n0_o + i, j, k);
		// 						else
		// 							result += norm(WF::psi[readInd1 + k])
		// 							* op.template operator() < R, opt > (j, n0_o + i, k);
		// 					}
		// 				}
		// 			}
		// 		}
		// 	}
		// }
		// return result * vol<R>();
		// return 0;
	// }
// };



// 	template <MODE M, REP R, OPTIMS opts> inline void prepare() const {}

// 	template<REP R, OPTIMS opt, typename Operator> inline double calc() const noexcept
// 	{
// 		constexpr Operator op;
// 	// auto op = COMPUTATION <double, false, Args...>::template getOperator<Operator>();
// 		double result = 0.0;
// 		// for (i = 0; i < n0_l; i++)
// 		// {
// 		// 	if constexpr (DIM == 1) result += norm(WF::psi[i]) * op. template operator() < R, opt > (i + n0_o);
// 		// 	else
// 		// 	{
// 		// 		readInd0 = i * n;
// 		// 		for (j = 0; j < n; j++)
// 		// 		{
// 		// 			readInd1 = readInd0 + j;
// 		// 			if constexpr (DIM == 2)
// 		// 			{
// 		// 				if constexpr (R == REP::X)
// 		// 					result += norm(WF::psi[readInd1])
// 		// 					* op. template operator() < R, opt > (i + n0_o, j);
// 		// 				else result += norm(WF::psi[readInd1])
// 		// 					* op. template operator() < R, opt > (j, i + n0_o);
// 		// 			}
// 		// 			else
// 		// 			{
// 		// 				readInd1 *= n;
// 		// 				for (k = 0; k < n; k++)
// 		// 				{
// 		// 					if constexpr (R == REP::X)
// 		// 						result += norm(WF::psi[readInd1 + k])
// 		// 						* op. template operator() < R, opt > (i + n0_o, j, k);
// 		// 					else result += norm(WF::psi[readInd1 + k])
// 		// 						* op. template operator() < R, opt > (j, i + n0_o, k);
// 		// 				}
// 		// 			}
// 		// 		}
// 		// 	}
// 		// }
// 		// return result * vol<R>();
// 		return 0;
// 	}
// };
