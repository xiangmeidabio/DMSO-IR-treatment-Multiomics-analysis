# Third-party Squidiff source

These files were copied from the local Squidiff analysis directory accompanying the manuscript. Squidiff is developed by Siyu He and collaborators. The upstream project is https://github.com/siyuh/Squidiff; the associated reproducibility repository is https://github.com/siyuh/Squidiff_reproducibility. Its source also acknowledges OpenAI guided-diffusion and Konpat's diffae code.

The local requirements file names Squidiff 1.0.8, but a matching upstream commit was not recorded. Accordingly, these files are identified by their source hashes, not described as an unmodified tagged upstream release. The local `train_squidiff.py` came from the manuscript's training folder. No neural-network algorithms were changed when preparing this repository.

The upstream local license file contained merge-conflict delimiters around the azizilab (2022) and Siyu He (2024) copyright notices. `LICENSE.upstream.txt` preserves that file exactly. `LICENSE` retains both notices and the BSD 3-Clause terms while removing only the conflict delimiters and duplicated final disclaimer sentence. These licenses apply to the vendored Squidiff code and do not impose a new license on the study's other source code.
