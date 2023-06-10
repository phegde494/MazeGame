# Maze Game

You can run the game via your favorite IDE. 

Just download the entire project repo to make sure you have access to the necessary JAR files.

This game uses Kruskal's minimum spanning tree algorithm to generate a random maze - with some adjustments, you can choose to add vertical or horizontal bias to the generation process. You can also view the maze being generated step by step via wall knockdowns if you so desire.

Generate as many random mazes as you want, or start over on the same maze.

You can use the arrow keys to manually move the player around and solve the maze yourself. The game will keep track of your incorrect moves. 

If you have trouble solving the maze, you can get the AI to solve it for you using either a BFS or DFS algorithm of your choice (DFS tends to work far better since the exit is always in the bottom right corner, far away from the start). 

You can also generate a heat map which displays each cell's distance from the start/end (your choice) in a gradient. To produce this gradient, I used a floodfill algorithm enabled via dynamic programming to update each cell's distance depending on its neighbors. 

As you or the AI moves across the grid, you can choose to display or hide all visited cells -> this visual can give you a good idea of the thinking process behind the search algorithm's decisions or help you keep track of your own progress.

To adjust the maze size, change the dimensions on line 926:

Currently, it says: `MazeWorld g = new MazeWorld(20, 10);` 

Adjust the values however you want up to a 100x60 maze.

Below is a demo of the game in which I demonstrate each of the features. Visited cells are highlighted in light blue, and their display can be turned on or off. First, I complete the maze manually. Then, I have the AI use DFS to solve the maze. Next, the AI uses BFS. Following this, I demonstrate the gradients and then create a new maze in four different ways: with a wall-knockdown animation, with horizontal bias, with vertical bias, and without any bias. You can follow the instructions on the right of the screen to do the same as you play.

**Have fun playing!**


https://github.com/phegde494/MazeGame/assets/48624928/4bbcd889-20f2-43b3-bc54-4a1c19cde9b0

