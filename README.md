This is a not-particularly-sophisticated CPU renderer that I wrote to cement my
understanding of the OpenGL pipeline's general structure, which it's ((very)
loosely) based on. It can draw a ~400-triangle sphere-with-holes with texturing
and basic Phong shading in 720x720 at about 30 FPS on my (fairly powerful)
laptop.

Working with some kind of GUI API to display the output would be out of scope
(the only dependencies are libc and libm), so the program streams video to
stdout as an undelimited sequence of raw RGBA images. If you have mpv
installed, you can use the `play.sh` wrapper script to watch the output in real
time.

The `LICENSE` file does not apply to `gear.pam`, `dna_vertices`, or
`dna_triangles`, which are not my original work.

