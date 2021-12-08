import matplotlib.pyplot as plt
   
problem_size = [128, 96, 640, 512, 1024, 2048, 4096]
cpu_time_efficient = [0.03125, 0.015625, 0.375, 0.46875, 1.78125, 8.1875, 39.3125]
memory_usage_efficient = [18879, 18185, 19642, 29412, 57932, 116468, 227652]
cpu_time_basic = [0.015625, 0.0, 0.25, 0.21875, 1.21875, 5.859375, 26.453125]
memory_usage_basic = [160761, 76591, 2489305, 2412718, 9555658, 38666124, 151490463]

# plt.plot(problem_size,cpu_time_efficient,label= 'efficient_memory')
# plt.plot(problem_size,cpu_time_basic,label= 'inefficient_memory')
# plt.title('Cpu_time Vs size')
# plt.xlabel('Problem_Size')
# plt.ylabel('Cpu_time')
# plt.legend()
# plt.show()

plt.plot(problem_size,memory_usage_efficient,label= 'efficient_memory')
plt.plot(problem_size,memory_usage_basic,label= 'inefficient_memory')
plt.title('Memory_usage Vs size')
plt.xlabel('Problem_Size')
plt.ylabel('Memory_Usage')
plt.legend()
plt.show()