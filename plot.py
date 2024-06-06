import matplotlib.pyplot as plt

# Step 1: Read Data from Text File
file_path = 'sum_output/sum_DYToLL.txt'
with open(file_path, 'r') as file:
    data = file.readlines()

# Step 2: Parse Data into X and Y Lists
x_values = []
y_values = []
for line in data:
    x, y = map(float, line.strip().split())  # Assuming columns are separated by whitespace
    x_values.append(x)
    y_values.append(y)

# Step 3: Plot the Data
plt.scatter(x_values, y_values)
plt.xlabel('X Axis Label')
plt.ylabel('Y Axis Label')
plt.title('Your Plot Title')
plt.grid(True)
plt.show()
