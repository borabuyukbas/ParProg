#include "plotting/plotter.h"
#include "io/image_parser.h"

#include <exception>

void Plotter::write_and_clear(){
    // create plot serial number string
    std::string serial_number_string = std::to_string(image_serial_number);
    while(serial_number_string.length() < 9){
        serial_number_string = "0" + serial_number_string;
    }

    std::string file_name = filename_prefix + "_" + serial_number_string + ".bmp";
    ImageParser::write_bitmap(output_folder_path / file_name, image);
    clear_image();
    image_serial_number += 1;
}

BitmapImage::BitmapPixel Plotter::get_pixel(std::uint32_t x, std::uint32_t y){
    return image.get_pixel(y, x);
}

void Plotter::highlight_position(Vector2d<double> position, std::uint8_t red, std::uint8_t green, std::uint8_t blue){
    // 
    if (!plot_bounding_box.contains(position)) return;

    // calculate pixel coordinates using relative width/height
    std::uint32_t x = ((position[0] - plot_bounding_box.x_min) / (plot_bounding_box.x_max - plot_bounding_box.x_min)) * plot_width;
    std::uint32_t y = ((position[1] - plot_bounding_box.y_min) / (plot_bounding_box.y_max - plot_bounding_box.y_min)) * plot_height;
    
    BitmapImage::BitmapPixel pixel(red, green, blue);

    //
    std::int32_t cross_width = 5;

    // draw cross
    for (int i = -cross_width; i <= cross_width; i++) {
        // set pixel on horizontal line
        image.set_pixel(y, x + i, pixel);
        // set pixel on vertical line
        image.set_pixel(y + i, x, pixel);
    }
}