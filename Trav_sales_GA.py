import numpy as np
import matplotlib.pyplot as plt
import random
import math

# Function to read city names and generate random coordinates
def load_cities(file_name):
    locations = []
    with open(file_name, 'r') as f:
        for line in f:
            city_name = line.strip()
            if city_name:
                # Generate random coordinates for cities
                longitude = random.uniform(60.0, 100.0)  
                latitude = random.uniform(8.0, 35.0)    
                locations.append((city_name, longitude, latitude))
            else:
                print(f"Skipping invalid line: {line}")
    return locations

# Function to calculate Euclidean distance between two cities
def city_distance(city_a, city_b):
    return math.sqrt((city_a[1] - city_b[1])**2 + (city_a[2] - city_b[2])**2)

# Function to calculate total route distance
def calculate_total_distance(path, locations):
    total_path_dist = 0
    for i in range(len(path) - 1):
        total_path_dist += city_distance(locations[path[i]], locations[path[i+1]])
        # Returning to start
    total_path_dist += city_distance(locations[path[-1]], locations[path[0]])   
    return total_path_dist

# Genetic Algorithm for TSP
def optimize_route(locations, num_individuals, max_generations, mutation_chance):
    # Initialize population randomly
    population = [random.sample(range(len(locations)), len(locations)) for _ in range(num_individuals)]

    avg_fitness_history = []
    best_route_history = []  # Tracking the best distance for each generation
    for gen in range(max_generations):
        # Calculate fitness for each individual (1 / total distance)
        fitness_scores = [1 / calculate_total_distance(route, locations) for route in population]
        all_distances = [calculate_total_distance(route, locations) for route in population]

        # Tournament selection for choosing parents
        selected_parents = []
        for _ in range(num_individuals):
            tournament = [random.randint(0, num_individuals-1) for _ in range(3)]
            winner = np.argmax([fitness_scores[i] for i in tournament])
            selected_parents.append(population[tournament[winner]])

        # Crossover: Order Crossover
        offspring_population = []
        for _ in range(num_individuals):
            parent1, parent2 = random.sample(selected_parents, 2)
            child = parent1[:len(parent1)//2] + [city for city in parent2 if city not in parent1[:len(parent1)//2]]
            offspring_population.append(child)

        # Mutation: Swap two cities with some probability
        for i in range(num_individuals):
            if random.random() < mutation_chance:
                idx1, idx2 = random.sample(range(len(offspring_population[i])), 2)
                offspring_population[i][idx1], offspring_population[i][idx2] = offspring_population[i][idx2], offspring_population[i][idx1]

        # Replace the old population with the new offspring
        population = offspring_population

        # Track the best route and its distance for this generation
        best_route_index = np.argmax(fitness_scores)
        best_route = population[best_route_index]
        best_route_distance = all_distances[best_route_index]
        best_route_history.append(best_route_distance)

        # Visualize the current best route
        plt.clf()
        plt.scatter([locations[i][1] for i in range(len(locations))], [locations[i][2] for i in range(len(locations))])
        for i in range(len(best_route) - 1):
            plt.plot([locations[best_route[i]][1], locations[best_route[i+1]][1]], [locations[best_route[i]][2], locations[best_route[i+1]][2]], 'b-')
        plt.plot([locations[best_route[-1]][1], locations[best_route[0]][1]], [locations[best_route[-1]][2], locations[best_route[0]][2]], 'b-')
        for i in range(len(locations)):
            plt.annotate(locations[i][0], (locations[i][1], locations[i][2]))
        plt.title(f'Generation {gen+1}')
        plt.pause(0.01)

        # Track average fitness for plotting
        avg_fitness_history.append(np.mean(fitness_scores))

    # Plot the average fitness improvement over generations
    plt.clf()
    plt.plot(avg_fitness_history)
    plt.xlabel('Generation')
    plt.ylabel('Average Fitness (1 / Distance)')
    plt.title('Improvement in Average Fitness')
    plt.show()

    # Plot the best route distance improvement
    plt.clf()
    plt.plot(best_route_history)
    plt.xlabel('Generation')
    plt.ylabel('Best Route Distance')
    plt.title('Improvement in Best Route Distance')
    plt.show()

    # Print final best route and its distance
    final_best_route = population[np.argmax(fitness_scores)]
    final_best_distance = min(best_route_history)
    print("Final Optimized Route:")
    for i in final_best_route:
        print(locations[i][0], end=" -> ")
    print(locations[final_best_route[0]][0])  # Close the loop
    print(f"Final Optimized Distance: {final_best_distance:.2f}")

# Main function
def main():
    cities_data = load_cities('india_cities.txt')   
    if not cities_data:
        print("No valid cities were loaded.")
        return

    num_individuals = 100
    max_generations = 100
    mutation_chance = 0.01
    optimize_route(cities_data, num_individuals, max_generations, mutation_chance)

if __name__ == '__main__':
    main()
