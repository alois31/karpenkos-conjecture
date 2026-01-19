# SPDX-FileCopyrightText: Alois Wohlschlager <wohlschlager@math.lmu.de>
# SPDX-License-Identifier: CC0-1.0
{
  development ? true,
}:
let
  pkgs =
    import
      (builtins.fetchTarball {
        url = "https://releases.nixos.org/nixos/unstable/nixos-25.11pre805967.62b852f6c674/nixexprs.tar.xz";
        sha256 = "sha256-0VFAX85s2QcFpcNVy7J4yn7j8Cn2IZADd1wsbT+1TYs=";
      })
      {
        overlays = [
          (final: prev: {
            nix = final.lixPackageSets.latest.lix;
            singular = prev.singular.overrideAttrs (old: {
              patches = old.patches or [ ] ++ [
                (final.fetchpatch {
                  url = "https://github.com/Singular/Singular/commit/926ff5410741911f311f3cef9ddfb4fd7c789bc3.diff";
                  hash = "sha256-bdfW/H2g4Mlog3JPz5sDhuniRwsbbAl1haKupVBgJVU=";
                })
                (final.fetchpatch {
                  url = "https://github.com/Singular/Singular/commit/0e16e23693fa02fbdc5e5208489b18a15fe30318.diff";
                  hash = "sha256-4mGqAT+3aw+tC+Cp9yUoP6RLzC2tPBwLa5a74qbsGig=";
                })
                (final.fetchpatch {
                  url = "https://github.com/Singular/Singular/commit/ac5186c27ddaa91a04cdd7a294fb91bdfc0ba281.diff";
                  hash = "sha256-O9CwclLhPWZx6zVKOxyF210IJUbD9vpUPmDcngVjszQ=";
                })
              ];
            });
          })
        ];
      };
in
pkgs.mkShell {
  strictDeps = true;
  nativeBuildInputs =
    [
      pkgs.cargo
      pkgs.rustc
      pkgs.singular
    ]
    ++ pkgs.lib.optionals development [
      pkgs.clippy
      pkgs.nil
      pkgs.nixfmt-rfc-style
      pkgs.reuse
      pkgs.rust-analyzer
      pkgs.rustfmt
    ];
}
