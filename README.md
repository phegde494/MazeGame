# Maze Game

You can run the game via your favorite IDE. 

Just download the entire project repo to make sure you have access to the necessary JAR files.

This game uses Kruskal's minimum spanning tree algorithm to generate a random maze - with some adjustments, you can choose to add vertical or horizontal bias to the generation process. You can also view the maze being generated step by step via wall knockdowns if you so desire.

Generate as many random mazes as you want, or start over on the same maze.

You can use the arrow keys to manually move the player around and solve the maze yourself. The game will keep track of your incorrect moves. 

If you have trouble solving the maze, you can get the AI to solve it for you using either a BFS or DFS algorithm of your choice (DFS tends to work far better since the exit is always in the bottom right corner, far away from the start). 

You can also generate a heat map which displays each cell's distance from the start/end (your choice) in a gradient. To produce this gradient, I used a breadth first search algorithm combined with dynamic programming to update each cell's distance depending on its neighbors. 

As you or the AI moves across the grid, you can choose to display or hide all visited cells -> this visual can give you a good idea of the thinking process behind the search algorithm's decisions or help you keep track of your own progress.

Have fun playing!
