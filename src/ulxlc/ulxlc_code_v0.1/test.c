#include <stdio.h>
#include <stdlib.h>
#include <sqlite3.h>

FTI = first transient index

int main(){
			erass 1												erass 2								erass 3								erass 4		
			#_ulx		#_not_ulx	P_ulx	P_not_ulx			#_FTI		P_transient				#_FTI		P_transient			#_FTI		P_transient
	sys1	569			404			0.569	0.404				214			0.214					342			0.342				110				0.110
	sys2
	sys3
	
	
	
	
	def quick_sample():
				for cycle in range(8):
					if cycle == 0:
						is_ulx = np.random.random() < system_cycle_1_P_ulx
					else:
						is_transient = np.random.random() < system_P_transient[cycle]
						if is_transient:
							return system_transient
							
				return system_persistent, is_ulx



	for repeat in range(1000):
		N_transients = 0
		N_alive_persistent = 0
		N_dead_persistent = 0

		for sys_n in range(500):
			for cycle in range(8):
				if cycle == 0:
					is_ulx = np.random.random() < system_cycle_1_P_ulx
				else:
					is_transient = np.random.random() < system_P_transient[cycle]
					if is_transient:
						return system_transient
						
			return system_persistent, is_ulx
	
	
	
	
class population:
	def __init__(self, systems):
		self.systems = [ulx1, ulx2, ...]
		self.size = len(self.systems)
		
		self.N_new_systems 				   = np.array([0,0,0,0,0,0,0,0])
		self.N_old_system_become_transient = np.array([None,0,0,0,0,0,0,0])
		self.N_delta_obs_ulx 			   = N_new_systems - N_old_system_become_transient
		self.N_observed_ulxs			   = cumsum(N_delta_obs_ulx)
		self.N_transients				   = N_new_systems + N_old_system_become_transient
		self.N_cum_transients 			   = cumsum(N_transients)
		self.N_total_systems			   = cumsum(N_new_systems)
	
	
	def run_sim(self):
		for cycle in range(8):
			new_systems = 0
			old_system_become_transient = 0
			
			for s in self.systems:
				if cycle == 0:
					is_ulx = s.is_ulx_first_cycle()
					new_systems += is_ulx
				else:
					if s.is_transient(cycle):
						if s.is_ulx:
							new_systems += 1
						else:
							N_old_system_become_transient += 1
					
			self.N_new_systems[cycle] = new_systems
			if cycle>0:
				self.N_old_system_become_transient[cycle] = old_system_become_transient


class ulx:
	def __init__(self):
		self.is_ulx = None
		self.P_cycle_1_ulx = 0.68
		self.P_transient = [None, 0.33, 0.46, 0.58, 0.69, 0.78, 0.99]
		self.transient_cycle = None
		
	def is_ulx_first_cycle(self):
		self.is_ulx = np.random.random() < self.P_cycle_1_ulx
		return self.is_ulx
		
	def is_transient(self, cycle):
		self.is_transient = np.random.random() < self.P_transient[cycle]
		if self.is_transient:
			self.is_ulx = not self.is_ulx
			self.transient_cycle = cycle
		return self.is_transient
			
			
	
				

	what we see:																											
																																			
							new_systems		old_system_become_transient		delta_obs_ulx		observed_ulxs		#_transients		#cum_transients		total_systems	
	cycle 1					50				None							+50			 		50					None				None           		50 				
	cycle 2					20				7								+13			 		63					27					27             		70				
	cycle 3					18				9								+9			 		72					27					54             		88				
	cycle 4					15				4								+11   		 		83			 	  	19					75             		103				
	cycle 5					13				16								-3			 		80					16					91             		116				
	cycle 6					12				2								+10			 		90					14					105            		128				
	cycle 7					8				1								+7			 		97					9					114            		136				
	
	
	
	N_new_systems 				  = observed
	N_old_system_become_transient = observed
	N_delta_obs_ulx 			  = N_new_systems - N_old_system_become_transient
	N_observed_ulxs				  =	cumsum(N_delta_obs_ulx)
	N_transients				  = N_new_systems + N_old_system_become_transient
	N_cum_transients 			  = cumsum(N_transients)
	N_total_systems				  = cumsum(N_new_systems)
	
	
	
	
	
	
	
	
	
	total_systems: 3
	
	erass 1:						sys1		sys2		sys3	sys4				cycle_transients 	cumulative_transients		systems_visible		total_systems_observed
	is_ulx?							ulx			not_ulx		ulx		not_ulx				2					None						2					2

	erass 2:																	
	observed as transient?			no			no			yes		no					1					1							1					2

	erass 3:
	observed as transient?			no			yes			...		no					1					2							2					3
	
	...
	
	erass 8:
	observed as transient?			no			...			...		no					0					2												
		
	
	
	
	system
	[0,1,0,1,1,1]
	
	if P_sup < 1000:
		do_p_sup
	if P_wind M < 1000:
		do_p_wind
		
	if P_sup < 1000 && p_wind < 1000
		do p_sup_&_p_wind
	
	
    return 0;

}


