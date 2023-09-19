let
  pkgs = import <nixpkgs> { };
in
pkgs.mkShell {
  name = "haskell stack shell";
  buildInputs = with pkgs; [
    stack
    haskell-language-server
  ];
}
