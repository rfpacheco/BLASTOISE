"""
BLASTOISE Module: Visual Aesthetics and Display Utilities
=========================================================

This module provides utility functions for creating visual elements in the BLASTOISE
pipeline's console output. It handles the formatting and display of messages and
ASCII art to enhance the user experience and provide visual feedback during pipeline
execution.

The module contains two main functions:
1. `print_message_box`: Creates a bordered box around text messages for emphasis
2. `blastoise_art`: Displays ASCII art of the Pokémon Blastoise, serving as a
   visual signature for the pipeline

These functions are used throughout the BLASTOISE pipeline to provide clear visual
cues about the progress and status of the analysis, making the console output more
readable and user-friendly.

Author: R. Pacheco
"""

BOX_WIDTH_DEFAULT = 80


def print_message_box(message: str, width: int = BOX_WIDTH_DEFAULT) -> None:
    """
    Print a message centered inside an ASCII art box for visual emphasis.

    This function creates a rectangular box using ASCII characters with the
    specified message centered inside it. The box has borders made of '=' characters
    and vertical bars '|' at the sides. The function automatically adjusts the box
    width if the message is longer than the specified width.

    Parameters
    ----------
    message : str
        Text to display centered inside the box.
    width : int, optional
        Number of '=' characters in the top/bottom borders.
        If smaller than the message length, it is automatically enlarged.
        Default is BOX_WIDTH_DEFAULT (80).

    Returns
    -------
    None
        The function prints directly to standard output and doesn't return a value.

    Examples
    --------
    >>> print_message_box("Hello World")

    |================================================================================|
    |                                  Hello World                                   |
    |================================================================================|

    >>> print_message_box("Short message", width=20)

    |====================|
    |   Short message    |
    |====================|
    """
    # Ensure the box is wide enough to fit the message
    width = max(width, len(message))  # ensure the message fits

    # Create the horizontal border using '=' characters
    border = "=" * width

    # Print the box with the message centered inside
    print()  # leading blank line
    print(f"|{border}|")  # top border
    print(f"|{message:^{width}}|")  # centered message
    print(f"|{border}|")  # bottom border
    print()  # trailing blank line


def blastoise_art() -> None:
    """
    Print ASCII art of Blastoise as a visual signature for the pipeline.

    This function displays a decorative ASCII art representation of the Pokémon
    Blastoise, which serves as a visual signature for the BLASTOISE pipeline.
    The art is printed to standard output and is typically used at the end of
    pipeline execution to provide a visually appealing completion indicator.

    The ASCII art includes water-like patterns surrounding the Blastoise figure,
    creating a thematic connection to the pipeline's name and purpose.

    Returns
    -------
    None
        The function prints directly to standard output and doesn't return a value.
    """
    print("""
                       o O       o O       o O       o O       o O
                     o | | O   o | | O   o | | O   o | | O   o | | O
                   O | | | | O | | | | O | | | | O | | | | O | | | | O
                  O-oO | | o   O | | o   O | | o   O | | o   O | | oO-o
                 O---o O o       O o       O o       O o      O o  O---o
                O-----O                                           O-----o
                o-----O           ⣠⣴⣾⣶⣿⣿⣶⣶⣶⣿⡟⠀⠀⠀                  o-----O
                 o---O          ⣠⣼⣿⣿⡿⣋⣠⣿⣿⣿⣿⣿⡶⢶⣶⣤⣤⣀⣤⣶⣿⣗⡤⠶⢦⡀⠀        o---O
                  o-O          ⣤⣾⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡿⠀⠘⢿⣿⣿⣿⡿⢁⣾⡀⠀⢀⡷⠀         o-O
                   O           ⠻⣿⠟⠛⠛⠻⠟⠛⠋⢹⣿⣿⣿⢣⡆⠀⠈⠛⠛⠋⠀⢻⣿⣿⣶⠟⣻⣦          O  
                  o-O             ⣴⡿⠁   ⢸⣿⣿⢇⣿⠃⠀⣠⣤⣤⣤⣤⣀⠉⠙⢁⣴⡿⠁         o-O    
                 o---O       ⢀⣾⣿⣿ ⢠⠀⠀⠀⣰⣤⣿⣿⢋⣬⡄⢀⣾⣿⣿⣿⣿⣿⣿⣧⠀⣿⣯⠀⠀        o---O    
                O-----O  ⢀⣾⣿⣿⣿⣿   ⣠⠻⠿⠿⠿⠿⣛⣵⣿⣿⣧⢸⣿⣿⣿⣿⣿⣿⣿⣿⣄⣿⣿⡆⠀       O-----o    
                O-----O⢀⣾⣿⣿⣿⣿    ⣠⣿⡀⢸⣿⣿⣿⣿⣿⣿⣿⠿⠆⠻⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣷⠀       O-----o  
                O---⣿⣿⣿⣿⣿⣿⣿  ⢀⣀⣴⣾⣿⣿⡇⣬⣭⣭⣭⣭⣭⣶⣶⣿⣷⡄⢈⣻⣿⣿⣿⣿⣿⣿⣿⣿⣿⠀       O-----o  
                o-⣿⣿⣿⣿⣿      ⠰⢾⣿⣿⣿⣿⡇⣿⣿⣿⣿⣿⣿⣿⣿⣿⡟⢐⣛⡻⣿⣿⣿⣿⣿⣿⠻⣿⣿⠀       o-----O
                 o-⣿⣿         ⠁⠀⣠⣶⣿⣷⢸⣿⣿⣿⣿⣿⣿⣿⡿⠿⠛⠋⡵⠿⢿⣿⣿⣿⢟⣄⢹⡏⠀        o---O
                  o-O          ⣰⣿⣿⣿⣿⣆⢲⣶⣶⣶⣶⣶⣶⣶⣿⢇⣷⣾⣿⡇⣟⣯⣶⣿⣿⡾⠀⠀         o-O  
                   O           ⣿⣿⣿⣿⣿⣿⣦⠹⣿⣿⣿⣿⣿⣿⣿⡜⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡾⠋          O
                  O-o          ⠘⣿⣿⣿⣿⣿⣿⣷⣬⠉⠿⣛⣻⣿⣯⣥⣹⣿⣿⣿⣿⣿⣿⣿⣿⣿⠀⠀         O-O
                 O---o        ⣠⣶⣿⣿⣿⣿⣿⣿⠿⠿⠦⠀⠀⠀⠉⠉⠁⠀⠹⣿⣿⣿⣿⣿⣿⣿⡿⠀⠀        O---o
                O-----o                         ⠀⠹⠿⠛⠿⣿⠟⠛⠛⠀        O-----o
                o-----O                                           o-----O
                 o---O o O       o O       o O       o O       o O o---O
                  o-Oo | | O   o | | O   o | | O   o | | O   o | | Oo-O
                   O | | | | O | | | | O | | | | O | | | | O | | | | O
                     O | | o   O | | o   O | | o   O | | o   O | | o
                       O o       O o       O o       O o       O o     
""")
