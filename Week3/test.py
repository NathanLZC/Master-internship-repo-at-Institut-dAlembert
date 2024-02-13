# define the constraints
# 平均压力约束??
def constraint_mean_pressure(p):
    return np.sum(p) / S - p_bar

# pressure bounds
bounds = [(0, None) for _ in range(len(p))]

# 初始压力分布假设??
p_initial = np.full((S,), p_bar)

# 定义优化问题??
def objective(p):
    x, y = np.meshgrid(np.linspace(-L/2, L/2, int(np.sqrt(len(p)))), np.linspace(-L/2, L/2, int(np.sqrt(len(p)))))
    x = x.flatten()
    y = y.flatten()
    return complementary_energy(p, x, y)

# 添加约束
constraints = [{'type': 'eq', 'fun': constraint_mean_pressure}]

# 进行优化
result = minimize(objective, p_initial, method='SLSQP', bounds=bounds, constraints=constraints)

if result.success:
    optimized_p = result.x
    print("Optimization successful.")
else:
    print("Optimization failed.")