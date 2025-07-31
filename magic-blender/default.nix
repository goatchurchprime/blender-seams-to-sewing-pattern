{ symlinkJoin, makeWrapper, lib, python3, blender, pkgs }:
let
  oldFreecad = (builtins.getFlake "github:nixos/nixpkgs/12a55407652e04dcf2309436eb06fef0d3713ef3").legacyPackages.${pkgs.hostPlatform.system}.freecad;
in
symlinkJoin {
  name = "blender-with-stuff";
  paths = [ blender ];
  buildInputs = [ makeWrapper ];
  postBuild = ''
    ls -lah $out/bin
    wrapProgram "$out/bin/blender" --suffix PYTHONPATH : ${
      lib.makeLibraryPath [ oldFreecad python3.pkgs.numpy ]
    }
  '';
}
