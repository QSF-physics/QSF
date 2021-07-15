
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
		// 	for (i = 0; i < local_m; i++) result += norm(WF::psi[i]);
		// }
		// else
		// {
		// 	for (ind i = 0; i < local_n; i++)
		// 	{
		// 		if constexpr (DIM == 1)
		// 			result += norm(WF::psi[i])
		// 			* op.template operator() < R, opt > (local_start + i);

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
		// 						* op.template operator() < R, opt > (local_start + i, j);
		// 					else
		// 						result += norm(WF::psi[readInd1])
		// 						* op.template operator() < R, opt > (j, local_start + i);
		// 				}
		// 				else
		// 				{
		// 					readInd1 = readInd1 * n;
		// 					for (ind k = 0; k < n; k++)
		// 					{
		// 						if constexpr (R == REP::X)
		// 							result += norm(WF::psi[readInd1 + k])
		// 							* op.template operator() < R, opt > (local_start + i, j, k);
		// 						else
		// 							result += norm(WF::psi[readInd1 + k])
		// 							* op.template operator() < R, opt > (j, local_start + i, k);
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



// template <typename... Args> struct CO_DIR_AVG : COMPUTATION <double, Args...>, _RBUFFER
// {
// 	static constexpr std::string_view name = "AVG_";
// 	// using type = double;
// 	// constexpr CO_DIR_AVG() :COMPUTATION<double, Args...>((Args::rep & ...)) {}

// 	template <REP R, typename WHEN>
// 	static constexpr bool canRun = COMPUTATION<double, Args...>::template goodRep<R>;

// 	template <MODE M, REP R, OPTIMS opts> inline void prepare() const {}

// 	template<REP R, OPTIMS opt, typename Operator> inline double calc() const noexcept
// 	{
// 		constexpr Operator op;
// 	// auto op = COMPUTATION <double, false, Args...>::template getOperator<Operator>();
// 		double result = 0.0;
// 		// for (i = 0; i < local_n; i++)
// 		// {
// 		// 	if constexpr (DIM == 1) result += norm(WF::psi[i]) * op. template operator() < R, opt > (i + local_start);
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
// 		// 					* op. template operator() < R, opt > (i + local_start, j);
// 		// 				else result += norm(WF::psi[readInd1])
// 		// 					* op. template operator() < R, opt > (j, i + local_start);
// 		// 			}
// 		// 			else
// 		// 			{
// 		// 				readInd1 *= n;
// 		// 				for (k = 0; k < n; k++)
// 		// 				{
// 		// 					if constexpr (R == REP::X)
// 		// 						result += norm(WF::psi[readInd1 + k])
// 		// 						* op. template operator() < R, opt > (i + local_start, j, k);
// 		// 					else result += norm(WF::psi[readInd1 + k])
// 		// 						* op. template operator() < R, opt > (j, i + local_start, k);
// 		// 				}
// 		// 			}
// 		// 		}
// 		// 	}
// 		// }
// 		// return result * vol<R>();
// 		return 0;
// 	}
// };
