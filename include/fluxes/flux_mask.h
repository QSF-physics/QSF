borVec<DIM>* mask;

/** to use this recalc gs and ev with Im::n = 128 and n=128 (or less)
 * otherwise the output will be unreadable in text editor*/
void previewMask(borVec<DIM>* mask, ind len, string_view name)
{
	ind index;
	int dir;
	constexpr double s = 0.001;
	logInfo("Previewing mask: %s", name.data());
	__logMPI("Node: %d\n", MPI::pID);

#if DIM >2
	for (k = 0; k < n; k++) //k will be the slowest index
	{
	#endif
	#if DIM >1
		for (j = 0; j < n; j++)
		{
		#if DIM==2
			__logMPI("j: %3td |", j);
		#elif DIM==3
			__logMPI("k: %3td j: %3td |", k, j);
		#endif
		#endif
			for (i = 0; i < len; i++)
			{
				index = i;
			#if DIM >1
				index = index * n + j;
			#endif
			#if DIM >2
				index = index * n + k;
			#endif

				dir = 0;
				if (mask[index][0] < -s) dir += 1; //left 
				else if (mask[index][0] > s)  dir += 2; //right
			#if DIM > 1
				if (mask[index][1] > s) dir += 4; // down
				else if (mask[index][1] < -s) dir += 8; //up
			#if DIM > 2
				if (mask[index][2] < -s) dir += 16; // back
				else if (mask[index][2] > s) dir += 32;//forw
			#endif
			#endif

				switch (dir) {
				case 0: __logMPI("%s", " "); break;
				case 1: __logMPI("%s", "←"); break;
				case 2: __logMPI("%s", "→"); break;
				case 4: __logMPI("%s", "↓"); break;
				case 8: __logMPI("%s", "↑"); break;
				case 5: __logMPI("%s", "↙"); break;
				case 6: __logMPI("%s", "↘"); break;
				case 9: __logMPI("%s", "↖"); break;
				case 10: __logMPI("%s", "↗"); break;
				case 16: __logMPI("%s", "-"); break;
				case 32:__logMPI("%s", "+"); break;
				default: __logMPI("%d", dir % 10); break;
				}
			}
			__logMPI("|\n");
		#if DIM >1
}
		__logMPI("\n");
	#endif
	#if DIM >2
	}
#endif
}

