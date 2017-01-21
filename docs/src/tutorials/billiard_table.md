# Setting up your Billiard Table
A billiard table `bt` is a vector of `Obstacle`s: `bt::Vector{Obstacle}`. 
The abstract Type `Obstacle` is the supertype of all objects that a particle may collide with.

There are some premade functions that construct well-known billiards, like the periodic Sinai billiard.
You can find all of them at the Standard Billiards page.

To create a custom billiard, you start with an empty Vector:
    bt = Obstacle[]
and then you create your obstacles one by one and add them to it. All obstacles that are already defined in the package
can be found at the Obstacles page.

